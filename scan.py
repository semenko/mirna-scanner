#!/usr/bin/env python
#
# Copyright (c) 2010 Nick Semenkovich <semenko@alum.mit.edu>
#
# Developed for the Nagarajan lab, Washington University in St. Louis (WUSTL)
# http://nagarajanlab.wustl.edu/
#
# This software is released under the MIT License <http://www.opensource.org/licenses/mit-license.php>
#
""" Perform multithreaded RNAhybrid alignments of miRNA orthologs to CD-masked ortholog clusters. """

from __future__ import generators
import cx_Oracle
import datetime
import logging
import math
import os
import Queue
import re
import scan_data
import socket
import subprocess
import threading
import time
import yappi
from optparse import OptionParser, OptionGroup
from bx.align import maf
from bx import interval_index_file

### -------------------------------
### Config Options
### -------------------------------


# Relative (or full) path to RNAhybrid
RNAHYBRID_PATH = 'rnahybrid/src/RNAhybrid'

# Same for filter_rnahybrid
FILTER_PATH = 'filter/filter_rnahybrid.pl'
FILTER_SETTINGS = 'filter/cutoff25'

# Database Settings
ORACLE_SERVER = '192.168.2.18'
ORACLE_SERVER = 'feservertest.wustl.edu'
ORACLE_DB = 'CHIPDB'
ORACLE_USERNAME = 'mirtarget'
ORACLE_PWFILE = '.oracle_password'

# Write result settings
ORACLE_WRITE_SERVER = '192.168.2.18'
ORACLE_WRITE_SERVER = 'feservertest.wustl.edu'
ORACLE_WRITE_DB = 'CHIPDB'
ORACLE_WRITE_USERNAME = 'nsemenkovich'
ORACLE_WRITE_PWFILE = '.oracle_write_password'

# Where to store log files of a run
LOG_PATH = 'results/'

# Where to read alignments from (we expect this path + GENEID/GENEID.maf.ordered
ALIGN_PATH = '/export/RAID_Drive_A/backup/Li-Wei/mtap/tba/seq/'

def SQLGenerator(cursor, arraysize = 2000):
    """ Don't fetchall (given ram) or fetchone (given net), so fetchmany! """
    
    while True:
        results = cursor.fetchmany(arraysize)
        if not results:
            break
        for result in results:
            yield result


REVCOMP = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}

def revComplement(sequence):
    """ Return the reverse complement of a sequence. """
    return ''.join([REVCOMP[bp] for bp in sequence[::-1]])


### ---------------------------------------------
### Exceptions
### ---------------------------------------------
class RNAHybridError(Exception):
    """ Thrown when RNAhybrid dies. """
    pass

class FilterError(Exception):
    """ Thrown when filter_rnahybrid.pl dies. """
    pass

### ---------------------------------------------


class RNAHybridThread(threading.Thread):
    """ The worker thread for RNAhybrid. """

    def __init__(self, thread_num, consume, input_queue, output_queue, species_index, weight_matrix, mrna_to_seq, mrna_to_exons):
        self.thread_num = thread_num
        self.consume = consume
        self.__input_queue = input_queue
        self.__output_queue = output_queue
        self.species_index = species_index # For looking up a weight
        self.weight_matrix = weight_matrix # Weights, for scaling exp(delta_g)
        self.logger = logging.getLogger('thread_%d' % thread_num)
        self.mrna_to_seq = mrna_to_seq
        self.mrna_to_exons = mrna_to_exons
        threading.Thread.__init__(self) # Remember to init our overridden base class

    def run(self):
        """ The workhorse of our thread. This works one one miRNA cluster -> homolog_cluster pair. """
        while True:
            self.consume.wait() # Wait until the "consume" event is True
            try:
                task = self.__input_queue.get(True, 5) # Block for at most 5 seconds on queue, otherwise die.
                # Check to make sure the output queue isn't full.
                # This shouldn't be the case, as Network Speed >>>> RNAHybrid Output Speed
                while self.__output_queue.full():
                    self.logger.error("Output queue is full! (Strange. Database issues?) Sleeping.")
                    time.sleep(1)

                micro_rna_id, micro_rna_cluster, hge_homologeneid, homolog_cluster = task # Unpack task from queue
                self.logger.info("Working on microrna: %s / homologene: %d" % (micro_rna_id, hge_homologeneid))
                # This is temporary storage for RNAhybrid results.
                filtered_rnahybrid = {}

                for local_taxid, microrna_seq in micro_rna_cluster:
                    try:
                        # WARNING: homolog_cluster can include more than one item, as the primary key for mRNA is
                        #          both mrc_geneid AND mrc_transcript_no

                        for mrc_geneid, mrc_transcript_no, exon_seq in homolog_cluster[local_taxid]:
                            ## Results is a set of tuples (all data in 5'->3' direction):
                            ## mfe, pos_from_3prime, target_3prime, target_bind, mirna_bind, mirna_5prime
                            ## ['-24.8', '742', '  G         GUAG  A     C', '   UGUGCAGCC    GC ACCUU ', '   AUAUGUUGG    UG UGGAG ', 'UUG            A  A     U']

                            # This spawns an rnahybrid subprocess (via rnahybrid(exon_seq, microrna_seq))
                            try:
                                filtered_rnahybrid.setdefault(local_taxid, []).append([microrna_seq,
                                                                                       mrc_geneid,
                                                                                       mrc_transcript_no,
                                                                                       rnahybrid(exon_seq, microrna_seq)])
                            except RNAHybridError:
                                self.logger.critical("RNAHybridError thrown in microrna: %s / homologene: %d" % (micro_rna_id, hge_homologeneid))
                            except FilterError:
                                self.logger.critical("FilterError thrown in microrna: %s / homologene: %d" % (micro_rna_id, hge_homologeneid))

                            print "JUST RAN on %d" % mrc_geneid
                            print "TAXID IS: %d" % local_taxid
                            print "MIRNA: %s" % microrna_seq
                            print "RAW SEQ: %s" % exon_seq
                                                            

                    except KeyError:
                        # No matching microRNA -> mRNA species pair.
                        pass
                print "LEN IS %d" % len(filtered_rnahybrid)

                ## We're done calling RNAhybrid and filtering it. Now let's compute scores and put results in the output queue.
                
                # Store result tuples of:
                # (micro_rna_id, local_taxid, microrna_seq, mrc_geneid, mrc_transcript_no, exon_seq, filtered_rnahybrid, unweighted_score, weighted_score)
                results = []

                # Loop over the human RNAhybrid runs, and compute their
                # weighted (using the MAF alignment files) and unweighted scores
                #
                # 7951 = human local tax id (TODO: Pull from dict, and shift dict to global element.)

                try:
                    for human_alignment in filtered_rnahybrid[7951]:
                        # Unpack values
                        microrna_seq, mrc_geneid, mrc_transcript_no, rnahybrid_data = human_alignment
                        
                        unweighted_score = 0.0
                        # Compute the unweighted sum(exp(delta_g)) score
                        for item in rnahybrid_data:
                            unweighted_score += math.exp(float(item[0]))
                            print "item here:"
                            print item
                            print "done"
                            print item[1]
                            print "item 3:"
                            print item[3]
                            # Match word boundaries and mark their positions in tuples
                            print [match.span() for match in re.finditer("\w+", item[3])]
                            print "word"
                            # The [match] iterable marks word boundary positions in tuples.
                            # This basically takes "ATG AGAGA GAGAGAG" and marks the start/stop pos of bases in the string.
                            # Then, map_local_to_global converts those ranges to their global exon position, e.g. (0, 3) -> (20500, 20503)
                            print _map_local_to_global_pos([match.span() for match in re.finditer("\w+", item[3])],
                                                           mrc_geneid, mrc_transcript_no, self.mrna_to_seq, self.mrna_to_exons)
                            
                        weighted_score = 0.0
                        try:
                            with open(ALIGN_PATH + str(mrc_geneid) + '/' + str(mrc_geneid) + '.maf.ordered', 'rb') as in_handle:
                                reader = maf.Reader(in_handle)
                                while True:
                                    pos = reader.file.tell()
                                    data = reader.next()
                                    if data is None:
                                        break
                                    for c in data.components:
                                        print c
                                        pass
                            # We lookup weights using: index = sum(2**species_index)
                            # self.species_index[]
                            # self.weight_matrix[]
                            in_handle.close()
                            print "align in %d" % mrc_geneid
                            self.logger.debug("Align found for: %d" % mrc_geneid)
                        except IOError:
                            print "no align in %d" % mrc_geneid
                            self.logger.debug("No align for: %d" % mrc_geneid)
                        print [micro_rna_id, hge_homologeneid,
                               local_taxid, mrc_geneid, mrc_transcript_no,
                               weighted_score, unweighted_score,
                               str(rnahybrid_data)]
                        print "---------------------------------------"
                        results.append([micro_rna_id, hge_homologeneid,
                                        local_taxid, mrc_geneid, mrc_transcript_no,
                                        weighted_score, unweighted_score,
                                        str(rnahybrid_data)])
                except KeyError:
                    self.logger.warning("No human taxonomy in microrna: %s / homologene: %d. Work was wasted." % (micro_rna_id, hge_homologeneid))

                self.__output_queue.put(results, True) # Block until a free slot is available.
            except Queue.Empty:
                # Nothing in the queue, so we're either done, or something bad has happened
                # and the main thread can't fill the Queue fast enough.
                self.logger.warning("Nothing in work_queue: Thread dying.")
                break


def rnahybrid(utr_seq, mirna_query):
    """ Execute RNAhybrid.
    
    Returns:
      Results set from run of RNAhybrid.
    """
    # Set storing RNAhybrid output
    # Value = (mfe, p-value, pos_from_3prime, target_3prime, target_bind, mirna_bind, mirna_5prime)
    results = []

    # We may later wish to add a p-value cutoff here, to help filter.pl out. Not clear.
    # RNAhybrid outputs:
    #    mfe, p_value, pos_from_3prime, target_3prime, target_bind, mirna_bind, mirna_5prime = line.split(':')[4:]
    #    results.add((float(mfe), float(p_value), int(pos_from_3prime), target_3prime, target_bind, mirna_bind, mirna_5prime))

    # First, run RNAhybrid
    process1 = subprocess.Popen([RNAHYBRID_PATH, '-c', '-s', '3utr_human', utr_seq, mirna_query],
                                bufsize = 1,
                                shell = False,
                                stdout = subprocess.PIPE)
    
    # and then run filter, piping the output in between.
    # *** filter is /undocumented/ and impossible to interpret. Not entirely clear what it does :(
    process2 = subprocess.Popen([FILTER_PATH, FILTER_SETTINGS],
                                bufsize = 1,
                                shell = False,
                                stdin = process1.stdout,
                                stdout = subprocess.PIPE,
                                stderr = subprocess.PIPE)
    stdoutdata2, stderrdata2 = process2.communicate()

    for line in stdoutdata2.split('\n')[:-1]:
        # We may wish to unpack these values, and store first two as int().
        # Filter output is \t separated, and looks like (after splitting):
        # ['-24.8', '742', '  G         GUAG  A     C', '   UGUGCAGCC    GC ACCUU ', '   AUAUGUUGG    UG UGGAG ', 'UUG            A  A     U']

        output = tuple(line.split('\t')[3:])
        # Filter to bindings of >=6 bases
        if re.search("[AUTGCN]{6,}", output[3]):
            results.append(output)

    # Sort the list by position. This can be costly :(
    results.sort(lambda x, y: cmp(int(x[1]), int(y[1])))

    if process1.returncode != None:
        raise RNAHybridError('An error occurred in executing RNAhybrid.')

    if process2.returncode != 0 or len(stderrdata2) != 0:
        raise FilterError('An error occurred filtering the RNAhybrid output: %s' % stderrdata2)


    return results


def sanity_overlap_check(in_dict):
    """ Sanity check to ensure start/stop positions don't overlap. """
    log_func = logging.getLogger('sanity_overlap_check')
    log_func.info("Validating bounds.")
    for value in in_dict.itervalues():
        start_old = 0
        end_old = 0
        for item in value:
            assert((item[0] > start_old) and (item[1] > start_old))
            assert((item[0] > end_old) and (item[1] > end_old))
            start_old = item[0]
            end_old = item[1]
    log_func.info("Validated.")


### ---------------------------------------------
### The main() function. Called first.
### ---------------------------------------------    
def main():
    """ Main execution. """

    # Begin timing execution
    starttime = time.time()

    # Check that RNAhybrid exists, and is executable
    assert(os.path.exists(RNAHYBRID_PATH) and os.access(RNAHYBRID_PATH, os.X_OK))
    # Same for Perl filter_rnahybrid script
    assert(os.path.exists(FILTER_PATH) and os.access(FILTER_PATH, os.X_OK))
    assert(os.path.exists(FILTER_SETTINGS))
    
    assert(os.path.exists(ALIGN_PATH))
    assert(os.path.exists(LOG_PATH))

    usage = "usage: %prog [OPTIONS] outfile"
    parser = OptionParser(usage)

    # This ignores hyperthreading pseudo-cores, which is fine since we hose the ALU
    parser.add_option("-j", help="Threads. We parallelize the invocations of RNAhybrid. [Default: # of CPU cores]",
                      default=os.sysconf('SC_NPROCESSORS_ONLN'), action="store", type="int", dest="threads")
    parser.add_option("-p", help="Profile. Invokes the yappi python profiling engine. Will slow down execution.",
                      default=False, action="store_true", dest="profileMe")

    group = OptionGroup(parser, "Range Settings (optional)")
    group.add_option("--start-num", help="What number miRNA ortholog group to start scanning from (inclusive).",
                      default=-1, action="store", type="int", dest="startNum")
    group.add_option("--stop-num", help="What number miRNA ortholog group to STOP scanning at (exclusive).",
                      default=-1, action="store", type="int", dest="stopNum")
    parser.add_option_group(group)

    
    (options, args) = parser.parse_args()

    # Sanity check range inputs
    range_given = False
    if options.startNum >= 0 or options.stopNum >= 0:
        if not (options.startNum >= 0 and options.stopNum >= 0):
            parser.error("If you specifiy a start/stop, you must specify both ends of the range!")
        if options.startNum >= options.stopNum:
            parser.error("Invalid scan range.")
        if options.startNum == 1:
            print "WARNING: This is a Python range, where lists are zero indexed! Are you sure you mean '1'?"
        range_given = True

    if len(args) == 0 and not range_given:
        parser.error("You must specify an outfile if you don't provide range options.\n\nTry -h for help.")

    # Set logging output, as flags passed were valid.
    # We log DEBUG and higher to log file, and write INFO and higher to console.
    datestamp = datetime.datetime.now().strftime("%m%d-%H%M")
    logfile = LOG_PATH + datestamp + '.' + socket.gethostname().split('.')[0]
    logging.basicConfig(filename = logfile, filemode = 'w',
                        format = '%(asctime)s: %(name)-25s: %(levelname)-8s: %(message)s',
                        level = logging.DEBUG)
    
    # define a Handler which writes INFO messages or higher to the sys.stderr
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    # set a format which is simpler for console use
    formatter = logging.Formatter('%(name)-25s: %(levelname)-8s: %(message)s')
    console.setFormatter(formatter)
    logging.getLogger('').addHandler(console)

    log_main = logging.getLogger('main')
    
    log_main.info('Logging to %s' % logfile)
    log_main.info('Starting run on %s' % socket.gethostname().split('.')[0])
    log_main.info('miRNA range: [%d, %d]' % (options.startNum, options.stopNum))

    if options.profileMe:
        log_main.info('Profiling enabled.')
        yappi.start()

    # Connect to Database
    oracle_password = open(ORACLE_PWFILE, 'r').readline()
    dbase = cx_Oracle.connect(ORACLE_USERNAME, oracle_password, ORACLE_SERVER + ':1521/' + ORACLE_DB)
    log_main.info("Connected to [read] %s (USER: %s)" % (dbase.dsn, ORACLE_USERNAME))

    # Connect to writing database, too
    oracle_write_password = open(ORACLE_WRITE_PWFILE, 'r').readline()
    dbase_write = cx_Oracle.connect(ORACLE_WRITE_USERNAME, oracle_write_password,
                                    ORACLE_WRITE_SERVER + ':1521/' + ORACLE_WRITE_DB)
    log_main.info("Connected to [write] %s (USER: %s)" % (dbase_write.dsn, ORACLE_WRITE_USERNAME))


    # Define the species sets we're looking at
    
    # Global TAXIDs. Local <-> Global <-> Org is available in CAFEUSER.LOCALTAXID_TO_ORG
    globalTaxMap = {'Hs': 9606, # Human
                    'Cf': 9615, # Dog
                    'Rn': 10116, # Rat
                    'Gg': 9031, # Chicken
                    'Mm': 10094} # Mouse
    revGlobalTaxMap = dict((v,k) for k, v in globalTaxMap.iteritems())

    localTaxMap = {'Hs': 7951, # Human
                   'Cf': 7959, # Dog
                   'Rn': 8385, # Rat
                   'Gg': 7458, # Chicken
                   'Mm': 8364} # Mouse
    revLocalTaxMap = dict((v,k) for k, v in localTaxMap.iteritems())

    UCSCTaxMap = {'hg18': 'Hs', # Human
                  'canFam2': 'Cf', # Dog
                  'rn4': 'Rn', # Rat
                  'galGal3': 'Gg', # Chicken
                  'mm9': 'Mm'} # Mouse
    revUCSCTaxMap = dict((v,k) for k, v in UCSCTaxMap.iteritems())

    shortToMAF = {'Hs': 'Human',
                  'Cf': 'Dog',
                  'Rn': 'Rat',
                  'Gg': 'Chicken',
                  'Mm': 'Mouse'}
    revShortToMAF = dict((v,k) for k, v in shortToMAF.iteritems())
    
    assert(len(globalTaxMap) == len(localTaxMap))
    assert(len(localTaxMap) == len(UCSCTaxMap))

    # Species Index for determining weights later
    # We lookup weights using: index = sum(2**species_index)
    species_index = {'hg18': 0,
                     'panTro2': 1,
                     'ponAbe2': 2,
                     'mm9': 3,
                     'rn4': 4,
                     'rheMac2': 5,
                     'monDom4': 6,
                     'bosTau4': 7,
                     'canFam2': 8,
                     'equCab2': 9}

    # Import our huge weight "matrix", which is a long list
    weight_matrix = scan_data.human_weight

    ### ---------------------------------------------
    # First, we get the microRNA clusters.
    mirna_queries = _get_microrna_data(dbase, range_given, options.startNum, options.stopNum, localTaxMap, UCSCTaxMap)

    ### ---------------------------------------------
    # Then, we get everything necessary for mRNA data.
    # Check out these functions for a description of these data structures.
    homologene_to_mrna = _get_homologene_to_mrna(dbase, globalTaxMap)
    mrna_to_seq = _get_mrna_to_seq(dbase, globalTaxMap)

    mrna_to_exons = _get_mrna_to_exons(dbase, globalTaxMap)
    sanity_overlap_check(mrna_to_exons)

    # Coding sequence excludes 3' and 5'
    #mrna_to_cds = _get_mrna_to_cds(dbase, globalTaxMap)
    #sanity_overlap_check(mrna_to_cds)
    ### ---------------------------------------------


    ### We do this by getting entire gene positions from PAP.GENE_COORDINATES, then
    ### selecting only exons using PAP.MRNA_COORDINATES.
    ### ---------------------------------------------
    ### First, start looping to populate our work_queue for threads to work
    ### Then, when filled, keep topping it off, while we periodically poll result_queue
    ### ---------------------------------------------

    work_queue_size = 500
    low_water_mark = 400 # Size at which we refill the work_queue
    fillup_increment = 50 # and how much to fill it by
    critical_low_mark = 100 # Size at which we stall threads, as work queue is too small!
    
    result_queue_size = 500
    high_water_mark = 100 # Size at which we send off data to collector (or write to disk)
    drain_increment = 50 # how many results to send off/write at a given time
    critical_high_mark = 400 # Site at which we stall threads, as result_queue is too big!

    assert(low_water_mark < work_queue_size)
    assert(critical_low_mark < low_water_mark)
    assert(high_water_mark < result_queue_size)
    assert(high_water_mark < critical_high_mark)
    assert(fillup_increment < (work_queue_size - low_water_mark))
    assert(drain_increment < (result_queue_size - high_water_mark))
    _out_of_work = False

    work_queue = Queue.Queue(maxsize = work_queue_size)  # A threadsafe producer/consumer Queue
    result_queue = Queue.Queue(maxsize = result_queue_size) # Same, but for product of RNAhybrid

    _first_run = True

    _work_queue_generator = _get_item_for_work_queue(dbase, mirna_queries, homologene_to_mrna, mrna_to_seq, mrna_to_exons)

    lol = []
    for _ in range(100):
        lol.append(_work_queue_generator.next())

    # WARNING: Queue calls can and /will/ block! Be careful!
    while True:
        log_main.info("(work|result|threads) = %s, %s, %s " %
                     (str(work_queue.qsize()), str(result_queue.qsize()), threading.active_count()-1))
        time.sleep(2)

        while (work_queue.qsize() < low_water_mark) and (_out_of_work == False):
            log_main.debug("Filling work queue.")
            # We call a generator to keep our work_queue mostly full
            for _ in range(fillup_increment):
                try:
                    work_queue.put(lol.pop())
#                    work_queue.put(_work_queue_generator.next())
                except StopIteration:
                    log_main.debug("Out of work: Waiting for work_queue to be emptied by threads ...")
                    _out_of_work = True
                except IndexError:
                    log_main.debug("Out of work: Waiting for work_queue to be emptied by threads ...")
                    _out_of_work = True
                    break
                                                            

        while result_queue.qsize() > high_water_mark:
            log_main.debug("Draining result_queue.")
            
            # These queue checks are imperfect, but should never be triggered.
            if result_queue.qsize() > critical_high_mark:
                log_main.error("Result queue is too big! Something is wrong!")
                consume_event.clear() # Stall threads.
            if (work_queue.qsize() < critical_low_mark) and (_out_of_work == False):
                log_main.error("Work queue is too small! Something is wrong!")
                consume_event.clear() # Stall threads.

            # Insert up to drain_increment elements into the database
            _write_results_to_db(dbase_write, result_queue, drain_increment)


        if _first_run:
            # First time we're running, so spawn threads
            log_main.info("Spawning %s threads." % str(options.threads))
            consume_event = threading.Event()
            for i in range(options.threads):
                thread = RNAHybridThread(i+1, consume_event, work_queue, result_queue, species_index, weight_matrix, mrna_to_seq, mrna_to_exons)
                thread.daemon = True
                thread.start()
            _first_run = False
            consume_event.set() # Let's go to work!

        if _out_of_work:
            # First, wait for results_queue to be completely emptied
            while (threading.active_count() - 1) > 0:
                _write_results_to_db(dbase_write, result_queue, 10)
            # At this point, all our threads have died.
            time.sleep(1) # A second to catch up, just in case.
            while result_queue.qsize() > 0:
                _write_results_to_db(dbase_write, result_queue, 1)
            log_main.info("All done!")
            break

        consume_event.set() # Let's go to work!
            

    log_main.info("Execution took: %s secs." % (time.time()-starttime))

    # Print output of profiling, if it was enabled.
    if options.profileMe:
        log_main.debug('Profiling data:')
        stats = yappi.get_stats()
        for stat in stats:
            log_main.debug(stat)
        yappi.stop()


def _write_results_to_db(dbase, result_queue, drain_increment):
    """ Write results to dbase. """
    write_cursor = dbase.cursor()
    write_cursor.prepare("INSERT INTO NSEMENKOVICH.RESULTS VALUES (:1, :2, :3, :4, :5, :6, :7, :8)")
    write_cursor.setinputsizes(None, None, None, None, None, None, None, cx_Oracle.CLOB)

    for _ in range(drain_increment):
        try:
            write_cursor.executemany(None, result_queue.get(True, 1))
        except Queue.Empty:
            time.sleep(1)
    
#        micro_rna_id, hge_homologeneid, local_taxid, mrc_geneid, mrc_transcript_no, weighted_score, unweighted_score, rnahybrid_data = result_queue.get(True, 5)

    dbase.commit()
    write_cursor.close()

def _get_microrna_data(dbase, range_given, startNum, stopNum, localTaxMap, UCSCTaxMap):
    """ Retreive microRNA data

        Args:
           dbase: a database object
           range_given: bool, take subset of miRNA ortholog groups, for
                        distributed computing
           startNum: start of range (implies range_given = True)
           stopNum: end of range (implies range_given = True)
           
        Returns:
          (dict) mirna_queries:
             keys: MMO_MATUREMIRID
             vals: ((localtaxid, MMO_MATURESEQ), (localtaxid, MMO_MATURESEQ) ...) """

    log_func = logging.getLogger('get_microrna_data')
    species_limit = ''.join(["'" + x + "'," for x in UCSCTaxMap.iterkeys()])[:-1]
    
    # Get list of distinct orthologous groups.
    mirna_id_dbcursor = dbase.cursor()
    mirna_id_dbcursor.execute("SELECT DISTINCT(MMO_MATUREMIRID) FROM LCHANG.MIRBASE_MIR_ORTHOLOG "
                              "WHERE MMO_MATURESEQ IS NOT NULL "
                              "AND MMO_SPECIES IN (" + species_limit + ") "
                              "ORDER BY 1")
    # validated:
    # mmo_maturemirid in ('hsa-miR-124', 'hsa-miR-1', 'hsa-miR-373','hsa-miR-155', 'hsa-miR-30a','hsa-let-7b')
    mirna_mirid_queries = mirna_id_dbcursor.fetchall()
    mirna_id_dbcursor.close()

    # (The list slicing is because fetchall() returns tuples.
    mirna_mirid_queries = [x[0] for x in mirna_mirid_queries]

    # We do some list slicing to get a range [if specified] of queries
    # This is so on a cluster, each machine can run over a selected few miRNAs
    if range_given:
        mirna_mirid_queries = mirna_mirid_queries[startNum:stopNum]
    log_func.info("Targeting %s miRNAs orthologous clusters" % len(mirna_mirid_queries))
    
    # Now that we have the list of miRNA orthologous clusters, make a dict of tuples
    # corresponding to the miRNA values.
    mirna_queries = {}
    # This will look like:
    # mirna_queries[MIR_ID] = ((localtaxid, MMO_MATURESEQ), (localtaxid, MMO_MATURESEQ) ...)

    mirna_dbcursor = dbase.cursor()
    # This is an inelegant way to select only the MIRIDs we want, but Oracle's lack of a
    # limit statement means this second query makes things more flexible with MySQL environments.
    mirna_dbcursor.execute("SELECT MMO_MATUREMIRID, MMO_SPECIES, MMO_MATURESEQ "
                           "FROM LCHANG.MIRBASE_MIR_ORTHOLOG "
                           "WHERE MMO_MATURESEQ IS NOT NULL "
                           "AND MMO_SPECIES IN ('hg18', 'canFam2', 'rn4', 'galGal3', 'mm9')")
    for row in SQLGenerator(mirna_dbcursor):
        if row[0] in mirna_mirid_queries:
            try:
                # Sanity check that miRNAs are actual sequences, free of 'N'
                assert(re.match("^[AUTGC]*$", row[2], re.IGNORECASE))
                # We convert MMO_SPECIES (e.g. 'canfam2') to a local taxid
                mirna_queries.setdefault(row[0], set()).add((localTaxMap[UCSCTaxMap[row[1]]], row[2]))
            except AssertionError:
                log_func.warning("Ignoring miRNA with 'N' or invalid char: %s, %s [mirid, species]"
                                 % (row[0], row[1]))
    mirna_dbcursor.close()

    # Sanity check that should never fail (unless you've broken the SQL queries)
    assert(len(mirna_mirid_queries) == len(mirna_queries))
    log_func.info("miRNA data retreived")
    
    return mirna_queries

def _get_homologene_to_mrna(dbase, globalTaxMap):
    """ Returns: homologene_to_mrna (dict) """
    # homologene_to_mrna: (Dict)
    #   key: hge_homologeneid # TODO: Double-check that not all HGE_GENEIDs map to MRC_GENEIDs
    #                         # Homologene has 244,950 rows, while the select has only 67,520
    #   val: set((mrc_geneid, mrc_transcript_no), (mrc_geneid, mrc_transcript_no), ...)
    log_func = logging.getLogger('get_homologene_to_mrna')
    homologene_to_mrna = {}
    species_limit = ''.join(["'" + str(x) + "'," for x in globalTaxMap.itervalues()])[:-1]
    log_func.info("Populating: homologene_to_mrna")
    mrna_dbase = dbase.cursor()
    mrna_dbase.execute("SELECT DISTINCT HGE_HOMOLOGENEID, MRC_GENEID, MRC_TRANSCRIPT_NO "
                       "FROM PAP.HOMOLOGENE, PAP.MRNA_COORDINATES "
                       "WHERE PAP.HOMOLOGENE.HGE_GENEID = PAP.MRNA_COORDINATES.MRC_GENEID "
                       "AND HGE_TAXID IN (" + species_limit + ")")
    for row in SQLGenerator(mrna_dbase):
        homologene_to_mrna.setdefault(row[0], set()).add((row[1], row[2]))
    mrna_dbase.close()
    log_func.info("Populated.")
    return homologene_to_mrna

def _map_local_to_global_pos(pos_ranges, mrc_geneid, mrc_transcript_no, mrna_to_seq, mrna_to_exons):
    """ The position returned by RNAhybrid is incorrect from a whole-genome perspective, as
        RNAhybrid only looks at the exon sequence. This computes positions with respect to
        introns, for use in MAF file parsing.
    
    Input:
       pos_range: A list of tuples of (START, STOP) following python (inclusive/exclusive) bounds.
       
    Returns:
       A list of position(s): [(start, stop), (start, stop), ...]
    """
    # TODO: Ponder if gcs_complement is important here.
    gcs_chromosome, gcs_taxid, gcs_localtaxid, gcs_complement, gcs_start, gcs_stop = mrna_to_seq[mrc_geneid]
    # val: [(mrc_start, mrc_stop), (mrc_start, mrc_stop), ...]
    try:
        exon_list = mrna_to_exons[(mrc_geneid, mrc_transcript_no)]
    except KeyError:
        self.logger.error("No exon mapping found! Something is terribly wrong! mrc_geneid: %s" % str(mrc_geneid))

    # Example data:
    # [(1, 7), (8, 9), (11, 15), (16, 19), (22, 26)] = pos_ranges
    # [(2137585, 2137803), (2138591, 2138794), (2138974, 2139015) ] = exon_list
    #
    # We're trying to map (1, 7) to (2137586, 2137592), bridging gaps in exons.
    # For example, mapping: (1, 7) on exon (100, 105), (105, 110) should return:
    # (101, 105), (105, 107) !!! Remember, python is (inclusive, exclusive)

    # Man, this function sucks. There has to be a more elegant way to map these values.
    # Do NOT do destrutive things here -- The list contents are inherited pointers!
    translated_positions = []
    exon_list_pointer = 0
    exon_offset = 0

    for l_start, l_end in pos_ranges:
        # Store our start and end positions
        s_pos = -1
        e_pos = -1
        while True:
            #print "list | offset: %d | %d" % (exon_list_pointer, exon_offset)
            #print "l_st | l_end: %d | %d" % (l_start, l_end)
            #print "s_pos | e_pos: %d | %d\n" % (s_pos, e_pos)
            exon = exon_list[exon_list_pointer]
            #print "exon is: %s" % str(exon)
            exon_len = exon[1] - exon[0]
            # Do we have our first start position yet?
            if (s_pos == -1):
                if exon_len >= (l_start - exon_offset):
                    s_pos = exon[0] + (l_start - exon_offset)
                else:
                    exon_list_pointer += 1
                    exon_offset += exon_len
                    continue # Skip rest of loop, go back to "while"
            if (e_pos == -1):
                if exon_len >= (l_end - exon_offset):
                    e_pos = exon[0] + (l_end - exon_offset)
                else:
                    exon_list_pointer += 1
                    exon_offset += exon_len
                    continue # Skip rest of loop, go back to "while"
            translated_positions.append((s_pos, e_pos))
            break

    # print "exon list:"
    # print exon_list
    # print "pos_list:"
    # print pos_ranges
    # print "translated!"
    # print translated_positions

    # Sanity check that (at least) the total number of basepairs match.
    assert(sum([pos[1]-pos[0] for pos in pos_ranges]) == sum([pos[1]-pos[0] for pos in translated_positions]))

    # Grab start/stop position of 
    return translated_positions

def _get_mrna_to_seq(dbase, globalTaxMap):
    """ Returns: mrna_to_seq (dict) """
    # mrna_to_seq: (Dict)
    #   key: mrc_geneid  # TODO: Double-check that mrc_transcript_no is NOT needed here
    #   val: (gcs_chromosome, gcs_taxid, gcs_localtaxid, gcs_complement, gcs_start, gcs_stop)
    log_func = logging.getLogger('get_mrna_to_seq')
    mrna_to_seq = {}
    species_limit = ''.join(["'" + str(x) + "'," for x in globalTaxMap.itervalues()])[:-1]
    log_func.info("Populating: mrna_to_seq")
    mrna_dbase = dbase.cursor()
    mrna_dbase.execute("SELECT DISTINCT(MRC_GENEID), GCS_CHROMOSOME, GCS_TAXID, "
                       "GCS_LOCALTAXID, GCS_COMPLEMENT, GCS_START, GCS_STOP "
                       "FROM PAP.MRNA_COORDINATES, PAP.GENE_COORDINATES "
                       "WHERE PAP.MRNA_COORDINATES.MRC_GENEID = PAP.GENE_COORDINATES.GCS_GENEID "
                       "AND GCS_TAXID IN (" + species_limit + ")")
    for row in SQLGenerator(mrna_dbase):
        assert(row[5] < row[6]) # Start < Stop
        mrna_to_seq[row[0]] = tuple(row[1:])
    mrna_dbase.close()
    log_func.info("Populated.")
    return mrna_to_seq


def _get_mrna_to_exons(dbase, globalTaxMap):
    """ Returns: mrna_to_exons (dict) """
    # mrna_to_exons: (Dict)
    #   key: (mrc_geneid, mrc_transcript_no)
    #   val: [(mrc_start, mrc_stop), (mrc_start, mrc_stop), ...]
    log_func = logging.getLogger('get_mrna_to_exons')
    mrna_to_exons = {}
    species_limit = ''.join(["'" + str(x) + "'," for x in globalTaxMap.itervalues()])[:-1]
    log_func.info("Populating: mrna_to_exons")
    mrna_dbase = dbase.cursor()
    mrna_dbase.execute("SELECT MRC_GENEID, MRC_TRANSCRIPT_NO, MRC_START, MRC_STOP "
                       "FROM PAP.MRNA_COORDINATES "
                       "WHERE MRC_GENEID IN (SELECT DISTINCT GCS_GENEID FROM PAP.GENE_COORDINATES "
                       "WHERE GCS_TAXID IN (" + species_limit + ")) "
                       "ORDER BY MRC_START")
    for row in SQLGenerator(mrna_dbase):
        assert(row[2] < row[3] + 1) # Start < Stop + 1
        mrna_to_exons.setdefault((row[0], row[1]), []).append((row[2], row[3]))
    mrna_dbase.close()
    log_func.info("Populated.")
    return mrna_to_exons


def _get_mrna_to_cds(dbase, globalTaxMap):
    """ Returns: mrna_to_cds (dict) """
    # mrna_to_cds: (Dict)
    #   key: (mrna_geneid, mrc_transcript_no)
    #   val: [(cds_start, cds_stop), (cds_start, cds_stop), ...]
    log_func = logging.getLogger('get_mrna_to_cds')
    mrna_to_cds = {}
    species_limit = ''.join(["'" + str(x) + "'," for x in globalTaxMap.itervalues()])[:-1]
    log_func.info("Populating: mrna_to_cds")
    mrna_dbase = dbase.cursor()
    mrna_dbase.execute("SELECT DISTINCT CDS_GENEID, CDS_TRANSCRIPT_NO, CDS_START, CDS_STOP "
                       "FROM PAP.CDS_COORDINATES "
                       "WHERE CDS_GENEID IN (SELECT DISTINCT GCS_GENEID FROM PAP.GENE_COORDINATES "
                       "WHERE GCS_TAXID IN (" + species_limit + ")) "
                       "ORDER BY CDS_START")
    for row in SQLGenerator(mrna_dbase):
        assert(row[2] < row[3] + 1) # Start < Stop + 1
        mrna_to_cds.setdefault((row[0], row[1]), []).append((row[2], row[3]))
    mrna_dbase.close()
    log_func.info("Populated.")
    return mrna_to_cds


def _get_item_for_work_queue(dbase, mirna_queries, homologene_to_mrna, mrna_to_seq, mrna_to_exons):
    """ A generator that yields one item at a time for the work queue. """

    seq_db = dbase.cursor() # TODO: Bind parameters in execute()
    
    # We work on one microRNA cluster at a time, and loop over all sequence clusters
    for micro_rna_id, micro_rna_cluster in mirna_queries.iteritems():
        # Loop over every Homologene=>(mrna_geneid, mrna_transcript_no)
        for hge_homologeneid, _mrna_list in homologene_to_mrna.iteritems():
            # homolog_cluster will get added to the input queue
            # This is a dict, rather than a simpler structure, so the threads can lookup
            # matching localtaxids for RNAhybrid comparisons.
            #
            # Key: gcs_localtaxid
            # Value: (mrc_geneid, mrc_transcript_no, exon_seq)
            homolog_cluster = {}

            # For each mRNA, lookup sequence (this is transcript no. agnostic)
            for mrc_geneid, mrc_transcript_no in _mrna_list:
                # Get sequence data for the entire gene
                gcs_chromosome, gcs_taxid, gcs_localtaxid, gcs_complement, gcs_start, gcs_stop = mrna_to_seq[mrc_geneid]

                # We add one to the sequence length since gcs_stop is strangely /inclusive/
                #  i.e. To select BP #70, gcs_start = 70 AND gcs_stop = 70
                seq_length = gcs_stop - gcs_start + 1

                # Extract the raw sequence from PAP.GENOMIC_SEQUENCE, which contains:
                #  > ges_loadid, ges_chromosome, ges_taxid, ges_localtaxid, ges_sequence, ges_masked_sequence
                #
                # Note: Avoid DBMS_LOB. functions, which get quirky with large result sets.
                # cx_Oracle can also handle CLOB objects directly, e.g.: row[0].read(1, 100)
                seq_db.execute("SELECT SUBSTR(GES_SEQUENCE," + str(gcs_start) + ", " + str(seq_length) + ") "
                               "FROM PAP.GENOMIC_SEQUENCE WHERE GES_TAXID = '" + str(gcs_taxid) + "' AND "
                               "GES_CHROMOSOME = '" + str(gcs_chromosome) + "'")
                seq_clob = seq_db.fetchone()[0]
                assert(seq_db.rowcount == 1)

                # Convert CLOB to STR for some easier handling
                raw_seq = str(seq_clob)
                assert(len(raw_seq) == seq_length)
                # Make sure we got actual base pairs, not '-'
                assert(re.match("^[ATGCN]*$", raw_seq, re.IGNORECASE))

                # Check if strand is complement, and take reverse complement.
                if gcs_complement:
                    raw_seq = revComplement(raw_seq)

                # Do some list slicing to get only the exon sequence.
                exons = []
                for exon_start, exon_stop in mrna_to_exons[(mrc_geneid, mrc_transcript_no)]:
                    # We have to offset these by gcs_start
                    # We also add one to compensate for python being [inclusive:exclusive]
                    exons.append(raw_seq[exon_start-gcs_start:exon_stop-gcs_start+1])
                exon_seq = ''.join(exons)

                # WARNING: We can have /multiple/ entries here, as the "primary key" for mRNA is
                # geneid AND transcript_no!
                homolog_cluster.setdefault(gcs_localtaxid, set()).add((mrc_geneid, mrc_transcript_no, exon_seq))

            yield (micro_rna_id, micro_rna_cluster, hge_homologeneid, homolog_cluster)
            
    seq_db.close()


if __name__ == "__main__":
    main()
