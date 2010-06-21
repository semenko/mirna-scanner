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

def rnahybrid(use_cache, pval_cutoff, organism_flag, entrez_geneid, utr_seq, mirna_query):
    """ Execute RNAhybrid.
    
    Returns:
      Results dictionary either from cache or a de-novo run of RNAhybrid.
    """

    # We use pickle to cache our results dict, and store it in a directory.
    # Hierarchy:
    #  mirna_query/entrez_geneid/organism_flag/organism_flag.pval_cutoff.cksum(utr_seq).out

    cacheDir = mirna_query + '/' + str(entrez_geneid) + '/' + organism_flag + '/'
    cacheKey = str(pval_cutoff) + '-' + str(zlib.adler32(utr_seq))

    if use_cache:
        results = load_cache('rnahybrid', cacheDir, cacheKey)
        if results:
            print "Got it!"
            # TODO: Add checks here for hash collisions.
            return results
        else:
            print "no results in cache :("

    # Dictionary Storing RNAhybrid output
    # Key = (entrez_geneid, organism_flag, mirna_sequence)
    # Value = (mfe, p-value, pos_from_3prime, target_3prime, target_bind, mirna_bind, mirna_5prime)
    results = {}

    process = subprocess.Popen([RNAHYBRID_PATH, '-c', '-p', pval_cutoff, '-s', '3utr_human', utr_seq, mirna_query],
                               shell=False,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)
    stdoutdata, stderrdata = process.communicate()
    for line in stdoutdata.split('\n')[:-1]:
        # ['command_line', '1957', 'command_line', '9', '-21.6', '0.597445', '966', 'A       G U', ' ACCCCGG G ', ' UGGGGCC C ', '        G  ']
        # we do this annoying unpack so the ints aren't stored as string
        # may not be necessary, as inputs may become strings in about 10 seconds ...
        mfe, p_value, pos_from_3prime, target_3prime, target_bind, mirna_bind, mirna_5prime = line.split(':')[4:]
        results[(int(entrez_geneid), organism_flag, mirna_seq)] = (float(mfe), float(p_value), int(pos_from_3prime), target_3prime, target_bind, mirna_bind, mirna_5prime)

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
        serial = pickle.load(open(CACHEPATH + module + '/' + cachedir + cachekey + '.cache', 'rb'))
        print serial
    except IOError:
        print "No cached data exists for these settings."
        return False
    except pickle.UnpicklingError:
        print "The cache was corruped. Regenerating."
        return False
### ---------------------------------------------

def tba():
    """ Execute TBA. """

    pass

def main():
    """ Main execution. """

    usage = "usage: %prog [OPTIONS]"
    parser = OptionParser(usage)

    # General settings
    parser.add_option("-o", "--includeOrthologs", help="Look up orthologs via MIRBASE_MIR_ORTHOLOG "
                      "and compute their bindings to QUERY as well. [Default: True]",
                      default=False, action="store_false", dest="useOrthologs")

    parser.add_option("-m", help="Specify a miRNA query sequence. If not specified, run over all "
                      "miRNA in the MIRBASE_MIR_ORTHOLOG database.",
                      default=False, action="store", type="string", dest="mirnaQuery")
    
    parser.add_option("-j", help="Threads. We perform pseudo-multithreading by parallelizing the "
                      "invocations of TBA. This is basically a spawned process limit. [Default: 2]",
                      default=2, action="store", type="int", dest="threads")

    parser.add_option("-c", "--usecache", help="Try to use local cache to retreive prior RNAhybrid"
                     "results. [Default: True]",
                     default=True, action="store_false", dest="useCache")

    # RNAhybrid Options
    group = OptionGroup(parser, "RNAhybrid Settings (optional)")

    group.add_option("-p", "--pvalue", help="p-value cutoff for results",
                     default=0.0, action="store", type="float", dest="pValue")

    parser.add_option_group(group)

    
    (options, args) = parser.parse_args()
    #if not options:
    #    parser.error("Try -h for help.")

    # Check that miRNA, if provided, is AUGC
    if (options.mirnaQuery):
        assert(re.match("^[AUTGCautgc]*$", args[-1]))

    # Instantiate SQL connection
    conn = MySQLdb.connect(host = "localhost",
                           user = "root",
                           passwd = "KleinersLaws",
                           db = "nagarajan")
    cursor = conn.cursor()

    # Organisms
    # Hs = Human
    # Pt = Chimp
    # Cf = Dog
    # Rn = Rat
    # Gg = Chicken
    # Mm = Mouse

    ### ------------------------------------------------------
    ### First, run RNAhybrid over the mirna_target for the given organism, otherwise all organisms.
    ### ------------------------------------------------------

    # This is a dict mapping short organisms ('Hs','Cf') to long titles
    # in the MIRBASE_MIR_ORTHOLOG database.
    organisms = {}

    # Should we look up orthologs to QUERY via MIRBASE_MIR_ORTHOLOG and scan them as well?
    if options.useOrthologs:
        pass
    else:
        # Don't look up orthlogs. Just use provided species and sequence.
        organism_flag = 'Hs'
        mirna_target = 'ACTGACATTTTGGGTCACA' # Fake target for test purposes.
        cursor.execute("SELECT ENTREZ_GENEID, UTR_SEQ FROM "
                       "T_PRM_UTRS_MIRTARGET WHERE UTR_COORDINATES "
                       "REGEXP \"^[0-9]*\_[0-9]*$\" AND TRANSCRIPT_NO = 0 "
                       "AND ORGANISM = '" + organism_flag +"' LIMIT 100")
        
        while True:
            val = cursor.fetchone()
            if val == None:
                break
            results = rnahybrid(options.useCache, options.pValue, organism_flag, val[0], val[1], mirna_target)
            print len(results)

            
if __name__ == "__main__":
    main()
