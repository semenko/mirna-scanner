#!/usr/bin/env python
#
# Copyright (c) 2010 Nick Semenkovich <semenko@alum.mit.edu>
#
# Developed for the Nagarajan lab, Washington University in St. Louis (WUSTL)
# http://nagarajanlab.wustl.edu/
#
# This software is released under the MIT License <http://www.opensource.org/licenses/mit-license.php>
#
""" HACKJOB: Distribute scan.py jobs across cluster computers. """

import cx_Oracle
import datetime
import math
import os
import re
import subprocess
import sys
import time

# Most settings are in scan.py. They should be given CLI flags.

CLUSTER_START = 1 # Inclusive
CLUSTER_STOP = 32 # Exclusive

# Database Settings
ORACLE_SERVER = 'feservertest.wustl.edu'
ORACLE_DB = 'CHIPDB'
ORACLE_USERNAME = 'mirtarget'
ORACLE_PWFILE = '.oracle_password'


def main():
    """ Main execution. """

    # Begin timing execution
    starttime = time.time()

    # Connect to Database
    oracle_password = open(ORACLE_PWFILE, 'r').readline()
    dbase = cx_Oracle.connect(ORACLE_USERNAME, oracle_password, ORACLE_SERVER + ':1521/' + ORACLE_DB)
    print "Connected to %s\n" % dbase.dsn

    UCSCTaxMap = {'hg18': 'Hs', # Human
                  'canFam2': 'Cf', # Dog
                  'rn4': 'Rn', # Rat
                  'galGal3': 'Gg', # Chicken
                  'mm9': 'Mm'} # Mouse
    
    species_limit = ''.join(["'" + x + "'," for x in UCSCTaxMap.iterkeys()])[:-1]
    
    # Get list of distinct orthologous groups.
    mirna_id_dbcursor = dbase.cursor()
    mirna_id_dbcursor.execute("SELECT COUNT(DISTINCT(MMO_MATUREMIRID)) FROM LCHANG.MIRBASE_MIR_ORTHOLOG "
                              "WHERE MMO_MATURESEQ IS NOT NULL "
                              "AND MMO_SPECIES IN (" + species_limit + ") ")
    mirna_mirid_queries = mirna_id_dbcursor.fetchall()
    mirna_id_dbcursor.close()

    number_of_mirna =  mirna_mirid_queries[0][0]
    comp_count = CLUSTER_STOP - CLUSTER_START

    increment = int(math.ceil(float(number_of_mirna) / comp_count))
    print "Distributing %s miRNA over %s machines [%s per machine]." % (number_of_mirna, comp_count, increment)


    #if range_given:
    #    mirna_mirid_queries = mirna_mirid_queries[startNum:stopNum]


    print "Distribution took: %s secs." % (time.time()-starttime)

if __name__ == "__main__":
    main()
