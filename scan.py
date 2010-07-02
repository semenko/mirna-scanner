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
import os
import Queue
import re
import socket
import subprocess
import sys
import threading
import time
#import yappi
import zlib
from optparse import OptionParser, OptionGroup

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

# Where to store results
RESULTS_PATH = 'results/'


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

    def __init__(self, thread_num, consume, input_queue, output_queue):
        self.thread_num = thread_num
        self.consume = consume
        self.__input_queue = input_queue
        self.__output_queue = output_queue
        threading.Thread.__init__(self) # Remember to init our overridden base class

    def run(self):
        """ The workhorse of our thread. """
        while True:
            self.consume.wait() # Wait until the "consume" event is True
            try:
                task = self.__input_queue.get(True, 5) # Block for at most 5 seconds on queue, otherwise die.
                # Check to make sure the output queue isn't full.
                # This shouldn't be the case, as Network Speed >>>> RNAHybrid Output Speed
                while self.__output_queue.full():
                    print "Output queue is full! (Strange. Network issues?) Sleeping."
                    time.sleep(1)
                # Value = (mfe, p-value, pos_from_3prime, target_3prime, target_bind, mirna_bind, mirna_5prime)
                micro_rna_id, micro_rna_cluster, hge_homologeneid, homolog_cluster = task # Unpack task from queue
                for local_taxid, microrna_seq in micro_rna_cluster:
                    try:
                        # WARNING: homolog_cluster can include more than one item, as the primary key for mRNA is
                        #          both mrc_geneid AND mrc_transcript_no
                        for mrc_geneid, mrc_transcript_no, exon_seq in homolog_cluster[local_taxid]:
                            ## Results is a set of tuples (all data in 5'->3' direction):
                            ## mfe, pos_from_3prime, target_3prime, target_bind, mirna_bind, mirna_5prime
                            ## ['-24.8', '742', '  G         GUAG  A     C', '   UGUGCAGCC    GC ACCUU ', '   AUAUGUUGG    UG UGGAG ', 'UUG            A  A     U']
                            results = rnahybrid(exon_seq, microrna_seq)
                    except KeyError:
                        # No matching microRNA -> mRNA species pair.
                        pass

                self.__output_queue.put(results, True) # Block until a free slot is available.
            except Queue.Empty:
                # Nothing in the queue, so we're either done, or something bad has happened
                # and the main thread can't fill the Queue fast enough.
                print "Nothing in work_queue: Thread %s dying." % self.thread_num
                break
                            
                

def rnahybrid(utr_seq, mirna_query):
    """ Execute RNAhybrid.
    
    Returns:
      Results set from run of RNAhybrid.
    """
    # Set storing RNAhybrid output
    # Value = (mfe, p-value, pos_from_3prime, target_3prime, target_bind, mirna_bind, mirna_5prime)
    results = set()

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
        print line.split('\t')[3:]
        results.add(tuple(line.split('\t')[3:]))
    
    if process1.returncode != None:
        raise RNAHybridError('An error occurred in executing RNAhybrid.')

    if process2.returncode != 0 or len(stderrdata2) != 0:
        raise FilterError('An error occurred filtering the RNAhybrid output: %s' % stderrdata2)


    return results


def sanity_overlap_check(in_dict):
    """ Sanity check to ensure start/stop positions don't overlap. """
    print "Validating bounds."
    for value in in_dict.itervalues():
        start_old = 0
        end_old = 0
        for item in value:
            assert((item[0] > start_old) and (item[1] > start_old))
            assert((item[0] > end_old) and (item[1] > end_old))
            start_old = item[0]
            end_old = item[1]
    print "Validated.\n"


### ---------------------------------------------
### The main() function. Called first.
### ---------------------------------------------    
def main():
    """ Main execution. """
    #yappi.start()

    # Begin timing execution
    starttime = time.time()

    # Check that RNAhybrid exists, and is executable
    assert(os.path.exists(RNAHYBRID_PATH) and os.access(RNAHYBRID_PATH, os.X_OK))
    # Same for Perl filter_rnahybrid script
    assert(os.path.exists(FILTER_PATH) and os.access(FILTER_PATH, os.X_OK))
    assert(os.path.exists(FILTER_SETTINGS))

    usage = "usage: %prog [OPTIONS] outfile"
    parser = OptionParser(usage)

    # This ignores hyperthreading pseudo-cores, which is fine since we hose the ALU
    parser.add_option("-j", help="Threads. We parallelize the invocations of RNAhybrid. [Default: # of CPU cores]",
                      default=os.sysconf('SC_NPROCESSORS_ONLN'), action="store", type="int", dest="threads")

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


    # Make our results path, and make sure we can write to it.
    if len(args) == 1:
        results_file = RESULTS_PATH + args[0]
        state_file = results_file + '.state'
    else:
        datestamp = datetime.datetime.now().strftime("%m%d-%H%M")
        results_file = RESULTS_PATH + datestamp + '.' + socket.gethostname().split('.')[0]

    print "\nWriting output to: %s" % results_file
    # Make our results path, if necessary
    if not os.path.exists(RESULTS_PATH):
        os.makedirs(RESULTS_PATH)
    open(results_file, 'w').close()

    # Connect to Database
    oracle_password = open(ORACLE_PWFILE, 'r').readline()
    dbase = cx_Oracle.connect(ORACLE_USERNAME, oracle_password, ORACLE_SERVER + ':1521/' + ORACLE_DB)
    print "Connected to %s\n" % dbase.dsn


    # Define the species sets we're looking at
    
    # Global TAXIDs. Local <-> Global <-> Org is available in CAFEUSER.LOCALTAXID_TO_ORG
    globalTaxMap = {'Hs': 9606, # Human
                    'Cf': 9615, # Dog
                    'Rn': 10116, # Rat
                    'Gg': 9031, # Chicken
                    'Mm': 10094} # Mouse

    localTaxMap = {'Hs': 7951, # Human
                   'Cf': 7959, # Dog
                   'Rn': 8385, # Rat
                   'Gg': 7458, # Chicken
                   'Mm': 8364} # Mouse

    UCSCTaxMap = {'hg18': 'Hs', # Human
                  'canFam2': 'Cf', # Dog
                  'rn4': 'Rn', # Rat
                  'galGal3': 'Gg', # Chicken
                  'mm9': 'Mm'} # Mouse
    
    assert(len(globalTaxMap) == len(localTaxMap))
    assert(len(localTaxMap) == len(UCSCTaxMap))


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

    stall_interval = 1 # How long (secs) to stall threads for (hopefully not necessary!)
    
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
    print "Work  | Result | Threads"
    while True:
        print " %s \t  %s \t  %s " % (str(work_queue.qsize()), str(result_queue.qsize()), threading.active_count()-1)
        time.sleep(2)

        while (work_queue.qsize() < low_water_mark) and (_out_of_work == False):
            print "Filling work queue."
            # We call a generator to keep our work_queue mostly full
            for _ in range(fillup_increment):
                try:
                    work_queue.put(lol.pop())
#                    work_queue.put(_work_queue_generator.next())
                except StopIteration:
                    print "Out of work: Waiting for work_queue to be emptied by threads ..."
                    _out_of_work = True
                except IndexError:
                    print "Out of work: Waiting for work_queue to be emptied by threads ..."
                    _out_of_work = True
                    break
                                                            

        while result_queue.qsize() > high_water_mark:
            print "Draining result_queue."
            
            # These queue checks are imperfect, but should never be triggered.
            if result_queue.qsize() > critical_high_mark:
                print "WARN: Result queue is too big! Something is wrong!"
                consume_event.clear() # Stall threads.
            if (work_queue.qsize() < critical_low_mark) and (_out_of_work == False):
                print "WARN: Work queue is too small! Something is wrong!"
                consume_event.clear() # Stall threads.

            # We reopen/close this handle, as NFS is sketchy.
            fh = open(results_file, 'a')
            for _ in range(drain_increment):
                # Send some results off to the collector, or write them to disk.
                fh.write(str(result_queue.get()))
            fh.close()


        if _first_run:
            # First time we're running, so spawn threads
            print "Spawning %s threads." % str(options.threads)
            consume_event = threading.Event()
            for i in range(options.threads):
                thread = RNAHybridThread(i, consume_event, work_queue, result_queue)
                thread.daemon = True
                thread.start()
            _first_run = False
            consume_event.set() # Let's go to work!

        if _out_of_work:
            # First, wait for results_queue to be completely emptied
            fh = open(results_file, 'a')
            while (threading.active_count() - 1) > 0:
                try:
                    fh.write(str(result_queue.get(True, 5)))
                except Queue.Empty:
                    pass
            # At this point, all our threads have died.
            time.sleep(1) # A second to catch up, just in case.
            while result_queue.qsize() > 0:
                fh.write(str(result_queue.get()))
            print "All done!"
            fh.close()
            break

        consume_event.set() # Let's go to work!
            

    print "Execution took: %s secs." % (time.time()-starttime)

#    stats = yappi.get_stats()
#    for stat in stats: print stat
#    yappi.stop()



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
    print "Targeting %s miRNAs orthologous clusters" % len(mirna_mirid_queries)
    
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
                print "\tIgnoring miRNA with 'N' or invalid char: %s, %s [mirid, species]" % (row[0], row[1])
    mirna_dbcursor.close()

    # Sanity check that should never fail (unless you've broken the SQL queries)
    assert(len(mirna_mirid_queries) == len(mirna_queries))
    print "miRNA data retreived.\n"
    
    return mirna_queries

def _get_homologene_to_mrna(dbase, globalTaxMap):
    """ Returns: homologene_to_mrna (dict) """
    # homologene_to_mrna: (Dict)
    #   key: hge_homologeneid # TODO: Double-check that not all HGE_GENEIDs map to MRC_GENEIDs
    #                         # Homologene has 244,950 rows, while the select has only 67,520
    #   val: set((mrc_geneid, mrc_transcript_no), (mrc_geneid, mrc_transcript_no), ...)
    homologene_to_mrna = {}
    species_limit = ''.join(["'" + str(x) + "'," for x in globalTaxMap.itervalues()])[:-1]
    print "Populating: homologene_to_mrna"
    mrna_dbase = dbase.cursor()
    mrna_dbase.execute("SELECT DISTINCT HGE_HOMOLOGENEID, MRC_GENEID, MRC_TRANSCRIPT_NO "
                       "FROM PAP.HOMOLOGENE, PAP.MRNA_COORDINATES "
                       "WHERE PAP.HOMOLOGENE.HGE_GENEID = PAP.MRNA_COORDINATES.MRC_GENEID "
                       "AND HGE_TAXID IN (" + species_limit + ")")
    for row in SQLGenerator(mrna_dbase):
        homologene_to_mrna.setdefault(row[0], set()).add((row[1], row[2]))
    mrna_dbase.close()
    print "Populated.\n"
    return homologene_to_mrna


def _get_mrna_to_seq(dbase, globalTaxMap):
    """ Returns: mrna_to_seq (dict) """
    # mrna_to_seq: (Dict)
    #   key: mrc_geneid  # TODO: Double-check that mrc_transcript_no is NOT needed here
    #   val: (gcs_chromosome, gcs_taxid, gcs_localtaxid, gcs_complement, gcs_start, gcs_stop)
    mrna_to_seq = {}
    species_limit = ''.join(["'" + str(x) + "'," for x in globalTaxMap.itervalues()])[:-1]
    print "Populating: mrna_to_seq"
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
    print "Populated.\n"
    return mrna_to_seq


def _get_mrna_to_exons(dbase, globalTaxMap):
    """ Returns: mrna_to_exons (dict) """
    # mrna_to_exons: (Dict)
    #   key: (mrc_geneid, mrc_transcript_no)
    #   val: [(mrc_start, mrc_stop), (mrc_start, mrc_stop), ...]
    mrna_to_exons = {}
    species_limit = ''.join(["'" + str(x) + "'," for x in globalTaxMap.itervalues()])[:-1]
    print "Populating: mrna_to_exons"
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
    print "Populated.\n"
    return mrna_to_exons


def _get_mrna_to_cds(dbase, globalTaxMap):
    """ Returns: mrna_to_cds (dict) """
    # mrna_to_cds: (Dict)
    #   key: (mrna_geneid, mrc_transcript_no)
    #   val: [(cds_start, cds_stop), (cds_start, cds_stop), ...]
    mrna_to_cds = {}
    species_limit = ''.join(["'" + str(x) + "'," for x in globalTaxMap.itervalues()])[:-1]
    print "Populating: mrna_to_cds"
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
    print "Populated.\n"
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
