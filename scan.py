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
import os
import cPickle as pickle
import cx_Oracle
import marshal
import threading
import Queue
import re
import subprocess
import time
import zlib
from optparse import OptionParser, OptionGroup

### -------------------------------
### Config Options
### -------------------------------


# Relative (or full) path to RNAhybrid
RNAHYBRID_PATH = 'rnahybrid/src/RNAhybrid'

# Where to store cold-cache data from RNAhybrid. This will be large (many GB).
COLDCACHE = 'cache/'

# Where to temporarily store rapidly churned cache data -- this should be a ram drive.
# (Hint: Mount /dev/shm somewhere.)
HOTCACHE = '/dev/shm/hot/'

# Database Settings
ORACLE_SERVER = 'feservertest.wustl.edu'
ORACLE_DB = 'CHIPDB'
ORACLE_USERNAME = 'mirtarget'
ORACLE_PWFILE = '.oracle_password'




def SQLGenerator(cursor, arraysize = 1000):
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


class RNAHybridThread(threading.Thread):
    """ The worker thread for RNAhybrid. """

    def __init__(self, thread_num, input_queue, output_queue):
        self.thread_num = thread_num
        self.__input_queue = input_queue
        self.__output_queue = output_queue
        threading.Thread.__init__(self) # Remember to init our overridden base class

    def run(self):
        """ The workhorse of our thread. """
        while True:
            task = self.__input_queue.get()
            if task is None:
                # Nothing in the queue, so we're either done, or something bad has happened
                # and the main thread can't fill the Queue fast enough.
                print "Nothing in work queue: Thread %s dying." % self.thread_num
                break
            else:
                # Check to make sure the output queue isn't full.
                # This shouldn't be the case, as Network Speed >>>> RNAHybrid Output Speed
                while self.__output_queue.full():
                    print "Output queue is full! (Strange. Network issues?) Sleeping."
                    time.sleep(1)
                # DO SOME WORK HERE
                print "Doing work."
                # Value = (mfe, p-value, pos_from_3prime, target_3prime, target_bind, mirna_bind, mirna_5prime)
                time.sleep(2)
                ## Results contains:
                ## (success_flag, invals, outvals)
                self.__output_queue.put(results, True) # Block until a free slot is available.
                

def rnahybrid(nocache, species, entrez_geneid, utr_seq, mirna_query):
    """ Execute RNAhybrid.
    
    Returns:
      Results set either from cache or a de-novo run of RNAhybrid.
    """

    cacheDir = mirna_query + '/' + str(entrez_geneid) + '/' + species + '/'
    cacheKey = str(zlib.adler32(utr_seq) & 0xffffffff) # Mask makes this positive.

    if not nocache:
        results = load_cache('rnahybrid', cacheDir, cacheKey)
        if results:
            print "\t\tRetreived from cache: %s (miRNA query)" % mirna_query
            # TODO: Add some addtl. validation of cache data?
            return results
    print "\t\tDe-novo run on: %s (miRNA query)" % mirna_query


    # Set storing RNAhybrid output
    # Value = (mfe, p-value, pos_from_3prime, target_3prime, target_bind, mirna_bind, mirna_5prime)
    results = set()

    # There's no reason to do a p-value cutoff here. We'd just be throwing away data that we might use later.
    process = subprocess.Popen([RNAHYBRID_PATH, '-c', '-s', '3utr_human', utr_seq, mirna_query],
                               shell=False,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)
    stdoutdata, stderrdata = process.communicate()
    for line in stdoutdata.split('\n')[:-1]:
        # This annoying unpack is to correctly store ints as ints, etc., and not everything as a string.
        mfe, p_value, pos_from_3prime, target_3prime, target_bind, mirna_bind, mirna_5prime = line.split(':')[4:]
        results.add((float(mfe), float(p_value), int(pos_from_3prime), target_3prime, target_bind, mirna_bind, mirna_5prime))

    if process.returncode != 0 or len(stderrdata) != 0:
        raise RNAHybridError('An error occurred in executing RNAhybrid.')

    save_cache('rnahybrid', cacheDir, cacheKey, results)
    return results

### ---------------------------------------------
### Caching code. # TODO: Use marshal instead.
### ---------------------------------------------
def save_cache(module, cachedir, cachekey, item):
    """ Store something in a pickle cache. """
    
    path = COLDCACHE + module + '/' + cachedir
    # Make sure all the directories exist.
    if not os.path.exists(path):
        os.makedirs(path)

    output = open(path + cachekey + '.cache', 'wb', -1)
    pickle.dump(item, output)
    output.close()
    return True

def load_cache(module, cachedir, cachekey):
    """ Try to retrieve something from a pickle cache. """
    try:
        return pickle.load(open(COLDCACHE + module + '/' + cachedir + cachekey + '.cache', 'rb'))
    except IOError:
        # No cached file exists
        return False
    # If this raises pickle.UnpicklingError, we have an error in the file itself. That's weird.
    # We don't catch it, since it's something you should look at.


def sanity_overlap_check(in_dict):
    """ Sanity check to ensure start/stop positions don't overlap. """
    print "Validating bounds."
    for value in in_dict.itervalues():
        start_old = 0
        end_old = 1
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

    # Begin timing execution
    starttime = time.time()

    # Check that RNAhybrid exists, and is executable
    #assert(os.path.exists(RNAHYBRID_PATH) and os.access(RNAHYBRID_PATH, os.X_OK))

    usage = "usage: %prog [OPTIONS]"
    parser = OptionParser(usage)

    # This ignores hyperthreading pseudo-cores, which is fine since we hose the ALU
    parser.add_option("-j", help="Threads. We parallelize the invocations of RNAhybrid. [Default: # of CPU cores]",
                      default=os.sysconf('SC_NPROCESSORS_ONLN'), action="store", type="int", dest="threads")

    parser.add_option("-n", "--nocache", help="Don't use local cache to retreive prior RNAhybrid"
                     "results. [Default: False]",
                     default=False, action="store_true", dest="noCache")

    group = OptionGroup(parser, "Range Settings (optional)")
    parser.add_option("--start-num", help="What number miRNA ortholog group to start scanning from (inclusive).",
                      default=-1, action="store", type="int", dest="startNum")
    parser.add_option("--stop-num", help="What number miRNA ortholog group to STOP scanning at (exclusive).",
                      default=-1, action="store", type="int", dest="stopNum")
    parser.add_option_group(group)

    
    (options, args) = parser.parse_args()
    if len(args) == 0:
        parser.error("Try -h for help.")

    # Sanity check range inputs
    range_given = False
    if options.startNum >= 0 or options.stopNum >= 0:
        if not (options.startNum >= 0 and options.stopNum >= 0):
            parser.error("If you specifiy a start/stop, you must specify both ends of the range!")
        if options.startNum >= options.stopNum:
            parser.error("Invalid scan range.")
        range_given = True

    # Connect to Database
    oracle_password = open(ORACLE_PWFILE, 'r').readline()
    dbase = cx_Oracle.connect(ORACLE_USERNAME, oracle_password, ORACLE_SERVER + ':1521/' + ORACLE_DB)
    print "Connected to %s\n" % dbase.dsn


    ### ---------------------------------------------
    # First, we get the microRNA clusters.
    mirna_queries = _get_microrna_data(dbase, range_given, options.startNum, options.stopNum)

    ### ---------------------------------------------
    # Then, we get everything necessary for mRNA data.
    # Check out these functions for a description of these data structures.
    homologene_to_mrna = _get_homologene_to_mrna(dbase)
    mrna_to_seq = _get_mrna_to_seq(dbase)
    
    mrna_to_exons = _get_mrna_to_exons(dbase)
    sanity_overlap_check(mrna_to_exons)

    mrna_to_cds = _get_mrna_to_cds(dbase)
    sanity_overlap_check(mrna_to_cds)
    ### ---------------------------------------------


    ### We do this by getting entire gene positions from PAP.GENE_COORDINATES, then
    ### selecting only exons using PAP.MRNA_COORDINATES.
    ### ---------------------------------------------
    ### First, start looping to populate our work_queue for threads to work
    ### Then, when filled, keep topping it off, while we periodically poll results_queue
    ### ---------------------------------------------

    work_queue_size = 1000
    low_water_mark = 900 # Size at which we refill the work_queue
    fillup_increment = 50 # and how much to fill it by
    critical_low_mark = 500 # Size at which we stall threads, as work queue is too small!
    
    result_queue_size = 1000
    high_water_mark = 100 # Size at which we send off data to collector (or write to disk)
    send_increment = 50 # how many results to send off/write at a given time
    critical_high_mark = 500 # Site at which we stall threads, as result_queue is too big!

    stall_interval = 1 # How long (secs) to stall threads for (hopefully not necessary!)
    
    assert(low_water_mark < work_queue_size)
    assert(fillup_increment < (work_queue_size - low_water_mark))
    _out_of_work = False

    work_queue = Queue.Queue(maxsize = work_queue_size)  # A threadsafe producer/consumer Queue
    # Contents:
    # (micro_rna_id, micro_rna_clusters, homologene_id, homolog_cluster)
    # Reminder: micro_rna_clusters = ((local_taxid, microrna_seq), ...)

    results_queue = Queue.Queue(maxsize = result_queue_size) # Same, but for product of RNAhybrid

    # WARNING: Queue calls can and /will/ block! Be careful!
    while True:
        while (work_queue.qsize() < low_water_mark) and (_out_of_work == False):
            # We call a generator to keep our work_queue mostly full
            for _ in range(fill_up_increment):
                try:
                    work_queue.add(_get_item_for_work_queue())
                except NoMoreWork:
                    print "Out of work: Waiting for work_queue to be emptied by threads ..."
                    _out_of_work = True

        while result_queue.qsize() > high_water_mark:
            if result_queue.qsize() > critical_high_mark:
                print "WARN: Results queue is too big! Something is wrong!"
                _do_insert_stall(work_queue, stall_interval)
            # Send some results off to the collector, or write them to disk.
            results_queue.get()
            

    print "Execution took: %s secs." % (time.time()-starttime)



def _get_microrna_data(dbase, range_given, startNum, stopNum):
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

    # This is a dict mapping short species tags ('Hs','Cf') to TAXIDs.
    # Note: These are /local/ TAXIDs via localtaxid_to_org (in cafeuser?)
    localTaxMap = {'Hs': 7951, # Human
                   'Pt': 7943, # Chimp
                   'Cf': 7959, # Dog
                   'Rn': 8385, # Rat
                   'Gg': 7458, # Chicken
                   'Mm': 8364} # Mouse
    
    UCSCTaxMap = {'hg18': 'Hs', # Human
                  'panTro2': 'Pt', # Chimp
                  'canFam2': 'Cf', # Dog
                  'rn4': 'Rn', # Rat
                  'galGal3': 'Gg', # Chicken
                  'mm9': 'Mm'} # Mouse

    # Get list of distinct orthologous groups.
    mirna_id_dbcursor = dbase.cursor()
    mirna_id_dbcursor.execute("SELECT DISTINCT(MMO_MATUREMIRID) FROM LCHANG.MIRBASE_MIR_ORTHOLOG "
                              "WHERE MMO_MATURESEQ IS NOT NULL "
                              "AND MMO_SPECIES IN ('mm9', 'rn4', 'canFam2', 'hg18') "
                              "ORDER BY 1")
    # validated:
    # where mmo_maturemirid in ('hsa-miR-124', 'hsa-miR-1', 'hsa-miR-373','hsa-miR-155', 'hsa-miR-30a','hsa-let-7b')
    #  and mmo_species in ('mm9', 'rn4', 'canFam2', 'hg18')
    mirna_mirid_queries = mirna_id_dbcursor.fetchall()
    mirna_id_dbcursor.close()

    # (The list slicing is because fetchall() returns tuples.
    mirna_mirid_queries = [x[0] for x in mirna_mirid_queries]

    # We do some list slicing to get a range [if specified] of queries
    # This is so on a cluster, each machine can run over a selected few miRNAs
    if range_given == True:
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
                           "AND MMO_SPECIES IN ('mm9', 'rn4', 'canFam2', 'hg18')")
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

def _get_homologene_to_mrna(dbase):
    """ Returns: homologene_to_mrna (dict) """
    # homologene_to_mrna: (Dict)
    #   key: hge_homologeneid # TODO: Double-check that not all HGE_GENEIDs map to MRC_GENEIDs
    #                         # Homologene has 244,950 rows, while the select has only 67,520
    #   val: set((mrc_geneid, mrc_transcript_no), (mrc_geneid, mrc_transcript_no), ...)
    homologene_to_mrna = {}
    print "Populating: homologene_to_mrna"
    mrna_dbase = dbase.cursor()
    mrna_dbase.execute("SELECT DISTINCT HGE_HOMOLOGENEID, MRC_GENEID, MRC_TRANSCRIPT_NO "
                       "FROM PAP.HOMOLOGENE, PAP.MRNA_COORDINATES "
                       "WHERE PAP.HOMOLOGENE.HGE_GENEID = PAP.MRNA_COORDINATES.MRC_GENEID "
                       "AND ROWNUM < 100 ")
    for row in SQLGenerator(mrna_dbase):
        homologene_to_mrna.setdefault(row[0], set()).add((row[1], row[2]))
    mrna_dbase.close()
    print "Populated.\n"
    return homologene_to_mrna


def _get_mrna_to_seq(dbase):
    """ Returns: mrna_to_seq (dict) """
    # mrna_to_seq: (Dict)
    #   key: mrc_geneid  # TODO: Double-check that mrc_transcript_no is NOT needed here
    #   val: (gcs_chromosome, gcs_taxid, gcs_localtaxid, gcs_complement, gcs_start, gcs_stop)
    mrna_to_seq = {}
    print "Populating: mrna_to_seq"
    mrna_dbase = dbase.cursor()
    mrna_dbase.execute("SELECT DISTINCT(MRC_GENEID), GCS_CHROMOSOME, GCS_TAXID, "
                       "GCS_LOCALTAXID, GCS_COMPLEMENT, GCS_START, GCS_STOP "
                       "FROM PAP.MRNA_COORDINATES, PAP.GENE_COORDINATES "
                       "WHERE PAP.MRNA_COORDINATES.MRC_GENEID = PAP.GENE_COORDINATES.GCS_GENEID "
                       "AND ROWNUM < 100")
    for row in SQLGenerator(mrna_dbase):
        assert(row[5] < row[6]) # Start < Stop
        mrna_to_seq[row[0]] = tuple(row[1:])
    mrna_dbase.close()
    print "Populated.\n"
    return mrna_to_seq


def _get_mrna_to_exons(dbase):
    """ Returns: mrna_to_exons (dict) """
    # mrna_to_exons: (Dict)
    #   key: (mrc_geneid, mrc_transcript_no)
    #   val: [(mrc_start, mrc_stop), (mrc_start, mrc_stop), ...]
    mrna_to_exons = {}
    print "Populating: mrna_to_exons"
    mrna_dbase = dbase.cursor()
    mrna_dbase.execute("SELECT MRC_GENEID, MRC_TRANSCRIPT_NO, MRC_START, MRC_STOP "
                       "FROM PAP.MRNA_COORDINATES "
                       "ORDER BY MRC_START")
    for row in SQLGenerator(mrna_dbase):
        assert(row[2] < row[3] + 1) # Start < Stop + 1
        mrna_to_exons.setdefault((row[0], row[1]), []).append((row[2], row[3]))
    mrna_dbase.close()
    print "Populated.\n"
    return mrna_to_exons


def _get_mrna_to_cds(dbase):
    """ Returns: mrna_to_cds (dict) """
    # mrna_to_cds: (Dict)
    #   key: (mrna_geneid, mrc_transcript_no)
    #   val: [(cds_start, cds_stop), (cds_start, cds_stop), ...]
    mrna_to_cds = {}
    print "Populating: mrna_to_cds"
    mrna_dbase = dbase.cursor()
    mrna_dbase.execute("SELECT DISTINCT MRC_GENEID, MRC_TRANSCRIPT_NO, CDS_START, CDS_STOP "
                       "FROM PAP.MRNA_COORDINATES, PAP.CDS_COORDINATES "
                       "WHERE PAP.MRNA_COORDINATES.MRC_GENEID = PAP.CDS_COORDINATES.CDS_GENEID "
                       "AND PAP.MRNA_COORDINATES.MRC_TRANSCRIPT_NO = PAP.CDS_COORDINATES.CDS_TRANSCRIPT_NO "
                       "AND ROWNUM < 100 "
                       "ORDER BY CDS_START")
    for row in SQLGenerator(mrna_dbase):
        assert(row[2] < row[3] + 1) # Start < Stop + 1
        mrna_to_cds.setdefault((row[0], row[1]), []).append((row[2], row[3]))
    mrna_dbase.close()
    print "Populated.\n"
    return mrna_to_cds


def _get_item_for_work_queue(mirna_queries, homologene_to_mrna, mrna_to_seq, mrna_to_exons):
    """ A generator that yields one item at a time for the work queue. """
    
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

                # TODO: Bind params, and stop reopening this damn cursor.
                seq_db = dbase.cursor()
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
                # Make sure we got actual base pairs, not '-' or 'N'
                assert(re.match("^[ATGC]*$", raw_seq, re.IGNORECASE))

                seq_db.close()

                # Check if strand is complement, and take reverse complement.
                if gcs_complement == True:
                    raw_seq = revComplement(raw_seq)

                # Do some list slicing to get only the exon sequence.
                exons = []
                for exon_start, exon_stop in mrna_to_exons[(mrc_geneid, mrc_transcript_no)]:
                    # We have to offset these by gcs_start
                    # We also add one to compensate for python being [inclusive:exclusive]
                    exons.append(raw_seq[exon_start-gcs_start:exon_stop-gcs_start+1])
                exon_seq = ''.join(exons)

                homolog_cluster[gcs_localtaxid] = (mrc_geneid, mrc_transcript_no, exon_seq)

            yield (micro_rna_id, micro_rna_cluster, homolog_cluster)


if __name__ == "__main__":
    main()
