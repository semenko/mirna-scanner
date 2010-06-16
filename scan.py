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

# Database Settings

import MySQLdb
import os
import cPickle as pickle
import subprocess
import sys
import time
from optparse import OptionParser, OptionGroup


def rnahybrid():
    """ Execute RNAhybrid. """

    # TODO: Allow user to change flags to RNAhybrid
    # str = "-c -s 3utr_human -t target -q query"

    conn = MySQLdb.connect(host = "localhost",
                           user = "root",
                           passwd = "KleinersLaws",
                           db = "nagarajan")
    cursor = conn.cursor()

    cursor.execute("SELECT UTR_SEQ, ENTREZ_GENEID FROM T_PRM_UTRS_MIRTARGET WHERE UTR_COORDINATES REGEXP \"^[0-9]*\_[0-9]*$\" AND TRANSCRIPT_NO = 0 AND ORGANISM = 'Hs'")

    mirna_seq = 'TGGGGCCGC'
    
    # Dictionary Storing RNAhybrid output
    # Key = (entrez_geneid, mirna_sequence)
    # Value = (mfe, p-value, pos_from_3prime, target_3prime, target_bind, mirna_bind, mirna_5prime)
    results = {}
    
    while True:
        val = cursor.fetchone()
        if val == None:
            break
        process = subprocess.Popen([RNAHYBRID_PATH, '-c', '-p', '0.8', '-s', '3utr_human', val[0], mirna_seq],
                                   shell=False,
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)
        stdoutdata, stderrdata = process.communicate()
        for line in stdoutdata.split('\n')[:-1]:
            # ['command_line', '1957', 'command_line', '9', '-21.6', '0.597445', '966', 'A       G U', ' ACCCCGG G ', ' UGGGGCC C ', '        G  ']
            # we do this annoying unpack so the ints aren't stored as string
            # may not be necessary, as inputs may become strings in about 10 seconds ...
            mfe, p_value, pos_from_3prime, target_3prime, target_bind, mirna_bind, mirna_5prime = line.split(':')[4:]
            results[(int(val[1]), mirna_seq)] = (float(mfe), float(p_value), int(pos_from_3prime), target_3prime, target_bind, mirna_bind, mirna_5prime)

        if process.returncode != 0:
            raise Exception('An error occurred in executing RNAhybrid.')

    print (len(results))

def save_cache(module, serial, item):
    """ Store something in a pickle cache. """
    
    path = 'cache/' + module + '/'
    if not os.path.exists(path):
        os.makedirs(path)

    output = open(path + serial + '.cache', 'wb', -1)
    pickle.dump(item, output)
    output.close()
    return True

def load_cache(module, serial):
    """ Try to retrieve something from a pickle cache. """
    try:
        serial = pickle.load(open('cache/' + module + '/' + serial + '.cache', 'rb'))
        print serial
    except IOError:
        print "No cached data exists for these settings."
        return False
    except pickle.UnpicklingError:
        print "The cache was corruped. Regenerating."
        return False
#    except:
#        print "Unable to retrieve from cache."
#        return False

def tba():
    """ Execute TBA. """

    pass

def main():
    """ Main execution. """

    # TODO: Add some assertion error checking?

    usage = "usage: %prog [OPTIONS]"
    parser = OptionParser(usage)

    # General settings
    parser.add_option("-j", help="Threads. We perform pseudo-multithreading by parallelizing the "
                      "invocations of TBA. This is basically a spawned process limit. [Default: 2]",
                      default=2, action="store", type="int", dest="threads")

    # RNAhybrid Options
    group = OptionGroup(parser, "RNAhybrid Settings (optional)")

    group.add_option("-m", "--longstring", help="Helpstring",
                     default=False, action="store_true", dest="variableAbc")

    parser.add_option_group(group)

    
    (options, args) = parser.parse_args()
    if len(args) == 0:
        parser.error("Try -h for help.")
    load_cache('mod', 'serial')
    save_cache('mod','serial','abc123')
#    rnahybrid()

            
if __name__ == "__main__":
    main()
