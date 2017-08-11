#!/usr/bin/env python

import sys
import tempfile
import os
import gzip
import shutil

from leafcutter_cluster import *

def main(options,libl):
    
    sort_junctions(libl, options)
    merge_junctions(options)
    get_numers(options)


if __name__ == "__main__":


    from optparse import OptionParser

    parser = OptionParser()

    parser.add_option("-j", "--juncfiles", dest="juncfiles",
                  help="text file with all junction files to be processed")

    parser.add_option("-o", "--outprefix", dest="outprefix", default = 'leafcutter',
                  help="output prefix (default leafcutter)")

    parser.add_option("-q", "--quiet",
                  action="store_false", dest="verbose", default=True,
                  help="don't print status messages to stdout")

    parser.add_option("-r", "--rundir", dest="rundir", default='./',
                  help="write to directory (default ./)")
    
    parser.add_option("-l", "--maxintronlen", dest="maxintronlen", default = 100000,
                  help="maximum intron length in bp (default 100,000bp)")

    parser.add_option("-m", "--minclureads", dest="minclureads", default = 30,
                  help="minimum reads in a cluster (default 30 reads)")

    parser.add_option("-p", "--mincluratio", dest="mincluratio", default = 0.001,
                  help="minimum fraction of reads in a cluster that support a junction (default 0.001)")

    (options, args) = parser.parse_args()

    if options.juncfiles == None:
        sys.stderr.write("Error: no junction file provided...\n")
        exit(0)
    
    # Get the junction file list
    libl = []
    for junc in open(options.juncfiles):
        junc = junc.strip()
        try:
            open(junc)
        except: 
            sys.stderr.write("%s does not exist... check your junction files.\n"%junc)
            exit(0)
        libl.append(junc)

    main(options, libl)
