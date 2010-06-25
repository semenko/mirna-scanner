#!/usr/bin/python
#
# Copyright (c) 2010 Nick Semenkovich <semenko@alum.mit.edu>
#
# Developed for the Nagarajan lab, Washington University in St. Louis (WUSTL)
# http://nagarajanlab.wustl.edu/
#
# This software is released under the MIT License <http://www.opensource.org/licenses/mit-license.php>
#

""" Perform multithreaded RNAhybrid alignments of miRNA orthologs to CD-masked ortholog clusters. """

### ---------------------------------------------
### Paths
### ---------------------------------------------

RNAHYBRID_PATH = 'rnahybrid/src/RNAhybrid'

# Where to store pickled cache data from RNAhybrid. This will be large (many GB).
COLDCACHE = 'cache/'

# Where to temporarily store output from pair-wise alignments (from all_bz).
# This has lots of data churn, and should be a ram drive. (Hint: Mount /dev/shm somewhere.)
HOTCACHE = '/dev/shm/hot/'


### ---------------------------------------------
## Database Settings
### ---------------------------------------------

ORACLE_SERVER = 'feservertest.wustl.edu'
ORACLE_DB = 'CHIPDB'
ORACLE_USERNAME = 'mirtarget'
ORACLE_PWFILE = '.oracle_password'


### ---------------------------------------------
### Imports
### ---------------------------------------------

import MySQLdb
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
### Pickling (caching) code.
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


### ---------------------------------------------
### The main() function. Called first.
### ---------------------------------------------    
def main():
    """ Main execution. """

    # Begin timing execution
    starttime = time.time()

    usage = "usage: %prog [OPTIONS]"
    parser = OptionParser(usage)

    # General settings
    #parser.add_option("-s", "--oneSpecies", help="Specify a species list for a query, such as "
    #                  "Hs,Pt [Default: Use all species.]",
    #                  default=False, action="store", type="string", dest="speciesList")

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

    # This is a dict mapping short species tags ('Hs','Cf') to TAXIDs.
    # Note: These are /local/ TAXIDs via localtaxid_to_org (in cafeuser?)
    speciesMap = {'Hs': 7951, # Human
                  'Pt': 7943, # Chimp
                  'Cf': 7959, # Dog
                  'Rn': 8385, # Rat
                  'Gg': 7458, # Chicken
                  'Mm': 8364} # Mouse

    # Either one or all species
    #speciesList = speciesMap.keys()
    #if options.oneSpecies:
    #    speciesList = [options.oneSpecies]

    # Instantiate SQL connection
    #dbase = MySQLdb.connect(host = "localhost", user = "root", passwd = "KleinersLaws", db = "nagarajan")

    oracle_password = open(ORACLE_PWFILE, 'r').readline()
    dbase = cx_Oracle.connect(ORACLE_USERNAME, oracle_password, ORACLE_SERVER + ':1521/' + ORACLE_DB)
    print "Connected to %s" % dbase.dsn

    # Get list of distinct orthologous groups.
    mirna_id_dbcursor = dbase.cursor()
    mirna_id_dbcursor.execute("SELECT DISTINCT(MMO_MATUREMIRID) FROM LCHANG.MIRBASE_MIR_ORTHOLOG "
                              "WHERE MMO_MATURESEQ IS NOT NULL ORDER BY 1")
    # validated:
    # where mmo_maturemirid in ('hsa-miR-124', 'hsa-miR-1', 'hsa-miR-373','hsa-miR-155', 'hsa-miR-30a','hsa-let-7b')
    #  and mmo_species in ('mm9', 'rn4', 'canFam2', 'hg18')
    mirna_mirid_queries = mirna_id_dbcursor.fetchall()
    mirna_id_dbcursor.close()
    
    # (The list slicing is because fetchall() returns tuples.
    mirna_mirid_queries = [x[0] for x in mirna_mirid_queries]

    # We do some list slicing here to get a range [if specified] of queries
    # This is so on a cluster, each machine can run over a selected few miRNAs    
    if range_given == True:
        mirna_mirid_queries = mirna_mirid_queries[options.startNum:options.stopNum]

    print "Targeting %s miRNAs orthologous clusters" % len(mirna_mirid_queries)

    
    # Now that we have the list of miRNA orthologous clusters, make a dict of tuples
    # corresponding to the miRNA values.
    mirna_queries = {}
    # This will look like:
    # mirna_queries[MIR_ID] = ((MMO_SPECIES, MMO_MATURESEQ), (MMO_SPECIES, MMO_MATURESEQ) ...)

    mirna_dbcursor = dbase.cursor()
    # This is an inelegant way to select only the MIRIDs we want, but Oracle's lack of a
    # limit statement means this second query makes things more flexible with MySQL environments.
    mirna_dbcursor.execute("SELECT MMO_MATUREMIRID, MMO_SPECIES, MMO_MATURESEQ "
                           "FROM LCHANG.MIRBASE_MIR_ORTHOLOG "
                           "WHERE MMO_MATURESEQ IS NOT NULL ORDER BY 1, 2")
    for row in mirna_dbcursor.fetchall():
        if row[0] in mirna_mirid_queries:
            # Sanity check that miRNAs are actual sequences
            try:
                assert(re.match("^[AUTGCautgc]*$", row[2]))
            except AssertionError:
                print row[2]
            mirna_queries.setdefault(row[0], set()).add((row[1], row[2]))
    mirna_dbcursor.close()

    # Sanity check that should never fail (unless you've broken the SQL queries)
    assert(len(mirna_mirid_queries) == len(mirna_queries))
    print "miRNA data retreived."

    ### ---------------------------------------------
    ### At this point, we have all the necessary miRNA clusters, so
    ### it's time to build /mRNA/ clusters from the database.
    ###
    ### We do this by getting the coding positions from @@@@@@@@@@@@, then
    ### selecting those regions from @@@@@@@@@@@.
    ###
    ### NOTE: Make sure to take the reverse compliment of the reverse strands.
    ###
    ### We preserve the gene id, and other meta data, so that later on, we can
    ### see if a region was a coding region, UTR, etc.
    ### ---------------------------------------------
    
    mrna_dbase = dbase.cursor()
    mrna_dbase.execute("SELECT MRC_GENEID, MRC_TRANSCRIPT_NO, MRC_START, MRC_STOP "
                       "FROM PAP.MRNA_COORDINATES WHERE ROWNUM <= 1000 "
                       "ORDER BY MRC_GENEID, MRC_TRANSCRIPT_NO, MRC_START")
    # We loop instead of fetchall(), as this is a huge result.
    # We could store CDS_COORDINTES as well at this step w/ a SQL JOIN, and might later, if RAM is available.
    mrna_coords = {}
    # This will be a dict:
    # Keys: (MRC_GENEID, MRC_TRANSCRIPT_NO)
    # Vals: [(MRC_START, MRC_STOP), (MRC_START, MRC_STOP), ...]
    # Note that vals is pre-sorted by MRC_START by the SQL SELECT statement.
    while True:
        row = mrna_dbase.fetchone()
        if row == None:
            break
        mrna_coords.setdefault((row[0], row[1]), []).append((row[2], row[3]))
    mrna_dbase.close()

    kold = None
    print "Validating bounds."
    for k in mrna_coords.iterkeys():
        for item in mrna_coords[k]:
            if k == kold:
                assert(item[0] > r2old)
                assert(item[1] > r2old)
                assert(item[0] > r3old)
                assert(item[1] > r3old)
                assert(item[0] < item[1])
                r2old = item[0]
                r3old = item[1]
            else:
                kold = k
                r2old = 0
                r3old = 1
    print "Validated."

    quit()

#    mirna_db_cursor.execute("SELECT ENTREZ_GENEID, UTR_SEQ FROM "
#                            "T_PRM_UTRS_MIRTARGET WHERE UTR_COORDINATES "
#                            "REGEXP \"^[0-9]*\_[0-9]*$\" AND TRANSCRIPT_NO = 0 "
#                            "AND ORGANISM = '" + species +"' "
#                            "LIMIT " + speed_limit)

    # Value = (mfe, p-value, pos_from_3prime, target_3prime, target_bind, mirna_bind, mirna_5prime)
    
    input_queue = Queue.Queue(maxsize = 1000) # A threadsafe producer/consumer Queue.
    output_queue = Queue.Queue(maxsize = 1000) # Same, but for product of RNAhybrid


    ### ---------------------------------------------
    ### First, start looping to populate our input_queue for threads to work
    ### Then, when filled, keep topping it off, while we periodically poll output_queue
    ### ---------------------------------------------

    # We work on one mirna cluster at a time, and loop over all sequence clusters
    for mmo_maturemirid in mirna_queries.iterkeys():
        """ Loop over mirna_queries dict handing out clusters of miRNAs to threads. """

        _, query_species, _ = mirna # Unpack: id, species, seq

        ## DATABASE SELECT HERE
        work_left = True
        while True:
            """ This loop runs until we've exhausted all sequences for the given mmo_mirid cluster. """
            
            while work_left == True:
                """ input_queue filling loop """

                # Fill work queue with mirna_cluster:sequence_cluster pairs
                # We loop over all sequence clusters, until our queue is full

                # keep advancing in database record
                mirna_sql = mirna_db_cursor.fetchone()
                if mirna_sql == None:
                    # We must be out of input miRNA ortholog groups
                    mirna_db_cursor.close()
                    work_left = False # Stop trying to fill the input_queue
                    break

        # Populate 

        print "\n\tUTR target: \t %s (entrez id)\n" % t_prm_sql[0]

        miRNAcursor = dbase.cursor()
        miRNAcursor.execute("SELECT DOM_MATURESEQ FROM DIST_ORTHOLOG_MIRNA WHERE DOM_LOCAL_TAXID = "
                            + str(speciesMap[species]) + " LIMIT " + speed_limit)

        while True:
            rna_sql = miRNAcursor.fetchone()
            if rna_sql == None:
                miRNAcursor.close()
                break

            # Put in queue
            queue.put((t_prm_sql[0], t_prm_sql[1], rna_sql[0]))
            #rnahybrid(options.noCache, species, t_prm_sql[0], t_prm_sql[1], rna_sql[0])
                
    print "\nAll done!\n"

    print "Execution took: %s secs." % (time.time()-starttime)

            
if __name__ == "__main__":
    main()
