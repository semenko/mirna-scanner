#!/usr/bin/python
#
# Copyright (c) 2010 Nick Semenkovich <semenko@alum.mit.edu>
#
# This software is released under the MIT License <http://www.opensource.org/licenses/mit-license.php>
#

""" Invokes RNAhybrid, passes to TBA. """

### ---------------------------------------------
### Paths to binaries
### ---------------------------------------------

RNAHYBRID_PATH = 'rnahybrid/src/RNAhybrid'
TBA_PATH = 'tba'

# Where to store pickled cache data from RNAhybrid. This will be large (many GB).
COLDCACHE = 'cache/'

# Where to temporarily store output from pair-wise alignments (from all_bz).
# This has lots of data churn, and should be a ram drive. (Hint: Mount /dev/shm somewhere.)
HOTCACHE = '/dev/shm/hot/'


### ---------------------------------------------
## Database Settings
### ---------------------------------------------



### ---------------------------------------------
### Imports
### ---------------------------------------------

import MySQLdb
import os
import cPickle as pickle
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


### ---------------------------------------------
### RNAhybrid Execution Module
### ---------------------------------------------
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

    if process.returncode != 0:
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
### TBA Execution Code for Alignments
### ---------------------------------------------
def tba():
    """ Execute TBA for multiple alignments. """

    species_data = "((((human chimp) (rat mouse)) dog) chicken)" # Is this the correct hierarchy?

    process = subprocess.Popen([TBA_PATH, '"' + species_data + '"', 'file here', 'file here'],
                               shell=False,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)
    stdoutdata, stderrdata = process.communicate()
    for line in stdoutdata.split('\n')[:-1]:
        pass


### ---------------------------------------------
### The main() function. Always called first.
### ---------------------------------------------    
def main():
    """ Main execution. """

    # Begin timing execution
    starttime = time.time()

    usage = "usage: %prog [OPTIONS]"
    parser = OptionParser(usage)

    # General settings
    parser.add_option("-s", "--oneSpecies", help="Specify a single species to query on, such as Hs "
                      "or Pt. [Default: Use all species.]",
                      default=False, action="store", type="string", dest="oneSpecies")

    parser.add_option("-m", help="Specify a miRNA query sequence. If not specified, run over all "
                      "miRNA in the MIRBASE_MIR_ORTHOLOG database.",
                      default=False, action="store", type="string", dest="mirnaQuery")
    
    #parser.add_option("-j", help="Threads. We perform pseudo-multithreading by parallelizing the "
    #                  "invocations of TBA. This is basically a spawned process limit. [Default: 2]",
    #                  default=2, action="store", type="int", dest="threads")

    parser.add_option("-n", "--nocache", help="Don't use local cache to retreive prior RNAhybrid"
                     "results. [Default: False]",
                     default=False, action="store_true", dest="noCache")

    # RNAhybrid Options
    # group = OptionGroup(parser, "RNAhybrid Settings (optional)")
    # parser.add_option_group(group)

    
    (options, args) = parser.parse_args()
    #if not options:
    #    parser.error("Try -h for help.")

    # Check that miRNA, if provided, is AUGC
    if (options.mirnaQuery):
        assert(re.match("^[AUTGCautgc]*$", args[-1]))
    
    # This is a dict mapping short species tags ('Hs','Cf') to TAXIDs.
    # Note: These are /local/ TAXIDs via localtaxid_to_org (in cafeuser?)
    speciesMap = {'Hs': 7951, # Human
                  'Pt': 7943, # Chimp
                  'Cf': 7959, # Dog
                  'Rn': 8385, # Rat
                  'Gg': 7458, # Chicken
                  'Mm': 8364} # Mouse

    # Either one or all species
    speciesList = speciesMap.keys()
    if options.oneSpecies:
        speciesList = [options.oneSpecies]

    # Instantiate SQL connection
    conn = MySQLdb.connect(host = "localhost",
                           user = "root",
                           passwd = "KleinersLaws",
                           db = "nagarajan")

    ### ------------------------------------------------------
    ### First, run RNAhybrid over the mirna_target for the given species, otherwise all species.
    ### ------------------------------------------------------

    # Build a Dictionary Storing RNAhybrid output
    # Key = (entrez_geneid, species, mirna_query)
    # Value = (mfe, p-value, pos_from_3prime, target_3prime, target_bind, mirna_bind, mirna_5prime)
    # results = {}

    speed_limit = "10"

    for species in speciesList:
        print "\nRunning RNAHybrid.\n\tspecies: \t %s" % species
        
        cursor = conn.cursor()
        cursor.execute("SELECT ENTREZ_GENEID, UTR_SEQ FROM "
                       "T_PRM_UTRS_MIRTARGET WHERE UTR_COORDINATES "
                       "REGEXP \"^[0-9]*\_[0-9]*$\" AND TRANSCRIPT_NO = 0 "
                       "AND ORGANISM = '" + species +"' "
                       "LIMIT " + speed_limit)

        while True:
            t_prm_sql = cursor.fetchone()
            if t_prm_sql == None:
                cursor.close()
                break

            print "\n\tUTR target: \t %s (entrez id)\n" % t_prm_sql[0]

            miRNAcursor = conn.cursor()
            miRNAcursor.execute("SELECT DOM_MATURESEQ FROM DIST_ORTHOLOG_MIRNA WHERE DOM_LOCAL_TAXID = "
                                + str(speciesMap[species]) + " LIMIT " + speed_limit)

            while True:
                rna_sql = miRNAcursor.fetchone()
                if rna_sql == None:
                    miRNAcursor.close()
                    break

                # This will return something. We just discard it for now (But it will build the cache.)
                # It clearly cannot fit in RAM, so we need some producer/consumer system.
                rnahybrid(options.noCache, species, t_prm_sql[0], t_prm_sql[1], rna_sql[0])
                
    print "\nAll done!\n"

    print "Execution took: %s secs." % (time.time()-starttime)

            
if __name__ == "__main__":
    main()
