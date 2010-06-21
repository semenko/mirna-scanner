#!/usr/bin/python
#
# Copyright (c) 2010 Nick Semenkovich <semenko@alum.mit.edu>
#
# This software is released under the MIT License <http://www.opensource.org/licenses/mit-license.php>
#

""" Invokes RNAhybrid, passes to TBA. """


# Binary Paths (Unnecessary if the execuatbles are in your $PATH)

RNAHYBRID_PATH = 'rnahybrid/src/RNAhybrid'
TBA_PATH = 'tba'

# Where to store cache data.
# You could make this /tmp/myfolder/, or make it a map to /dev/shm for in-RAM caching.
CACHEPATH = 'cache/'

# Database Settings

import MySQLdb
import os
import cPickle as pickle
import re
import subprocess
import sys
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

def rnahybrid(nocache, species, entrez_geneid, utr_seq, mirna_query):
    """ Execute RNAhybrid.
    
    Returns:
      Results set either from cache or a de-novo run of RNAhybrid.
    """

    cacheDir = mirna_query + '/' + str(entrez_geneid) + '/' + species + '/'
    cacheKey = str(zlib.adler32(utr_seq))

    if not nocache:
        results = load_cache('rnahybrid', cacheDir, cacheKey)
        if results:
            print "\t\tRetreived from cache: %s (entrez id)" % entrez_geneid
            # TODO: Add checks here for hash collisions.
            return results
    print "\t\tDe-novo run on: %s (entrez id)" % entrez_geneid


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
        # ['command_line', '1957', 'command_line', '9', '-21.6', '0.597445', '966', 'A       G U', ' ACCCCGG G ', ' UGGGGCC C ', '        G  ']
        # we do this annoying unpack so the ints aren't stored as string
        # may not be necessary, as inputs may become strings in about 10 seconds ...
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
    
    path = CACHEPATH + module + '/' + cachedir
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
        cachevals = pickle.load(open(CACHEPATH + module + '/' + cachedir + cachekey + '.cache', 'rb'))
        return cachevals
    except IOError:
        # No cached file exists
        return False
    # If this raises pickle.UnpicklingError, we have an error in the file itself. That's weird.

### ---------------------------------------------

def tba():
    """ Execute TBA. """

    pass

def main():
    """ Main execution. """

    usage = "usage: %prog [OPTIONS]"
    parser = OptionParser(usage)

    # General settings
    parser.add_option("-s", "--oneSpecies", help="Specify a single species to query on, such as Hs "
                      "or Pt. [Default: Use all species.]",
                      default=False, action="store", type="string", dest="oneSpecies")

    parser.add_option("-m", help="Specify a miRNA query sequence. If not specified, run over all "
                      "miRNA in the MIRBASE_MIR_ORTHOLOG database.",
                      default=False, action="store", type="string", dest="mirnaQuery")
    
    parser.add_option("-j", help="Threads. We perform pseudo-multithreading by parallelizing the "
                      "invocations of TBA. This is basically a spawned process limit. [Default: 2]",
                      default=2, action="store", type="int", dest="threads")

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
    
    # This is a dict mapping short species tags ('Hs','Cf') to long titles
    # in the MIRBASE_MIR_ORTHOLOG database and TAXIDs.
    # Note: These are /local/ TAXIDs via localtaxid_to_org (in cafeuser?)
    speciesMap = {'Hs': ('', 7951), # Human
                  'Pt': ('', 7943), # Chimp
                  'Cf': ('', 7959), # Dog
                  'Rn': ('', 8385), # Rat
                  'Gg': ('', 7458), # Chicken
                  'Mm': ('', 8364)} # Mouse


    # Either one or all species
    speciesList = speciesMap.keys()
    if options.oneSpecies:
        speciesList = [options.oneSpecies]

    # Instantiate SQL connection
    conn = MySQLdb.connect(host = "localhost",
                           user = "root",
                           passwd = "KleinersLaws",
                           db = "nagarajan")
    cursor = conn.cursor()

    ### ------------------------------------------------------
    ### First, run RNAhybrid over the mirna_target for the given species, otherwise all species.
    ### ------------------------------------------------------
    mirna_query = 'ACTGACATTTTGGGTCACA' # Fake target for test purposes.

    # Build a Dictionary Storing RNAhybrid output
    # Key = (entrez_geneid, species, mirna_query)
    # Value = (mfe, p-value, pos_from_3prime, target_3prime, target_bind, mirna_bind, mirna_5prime)
    results = {}
                        
    for species in speciesList:
        print "\nRunning RNAHybrid.\n\tspecies: \t %s" % species
        print "\tmiRNA Query: \t %s\n" % mirna_query
    
        cursor.execute("SELECT ENTREZ_GENEID, UTR_SEQ FROM "
                       "T_PRM_UTRS_MIRTARGET WHERE UTR_COORDINATES "
                       "REGEXP \"^[0-9]*\_[0-9]*$\" AND TRANSCRIPT_NO = 0 "
                       "AND ORGANISM = '" + species +"'")

        while True:
            val = cursor.fetchone()
            if val == None:
                break
            # Full results storage will kill RAM. This needs to be a deque, or something.
            results[(val[0], species, mirna_query)] = rnahybrid(options.noCache, species, val[0], val[1], mirna_query)

    print "All done!\nResults is: %s" % len(results)

            
if __name__ == "__main__":
    main()
