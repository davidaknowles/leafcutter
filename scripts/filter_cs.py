#!/usr/bin/env python
"""
filter_cs.py

Filter alignments in a SAM file by the CIGAR string and mapping quality.
"""

import sys
import argparse
import fileinput

def parse_args():
    parser = argparse.ArgumentParser(
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("--min_intron_length",
            type=int, action="store", dest="min_intron_length", default=50,
            help="Alignments spanning introns shorter than this are ignored.")

    parser.add_argument("--min_overlap",
            type=int, action="store", dest="min_overlap", default=6,
            help="The 5' and 3' edges of a read that spans 2 exons"
                 " must both overlap this many bases into each exon.")

    parser.add_argument("--min_quality",
            type=int, action="store", dest="min_quality", default=10,
            help="Alignments with less than this quality are ignored.")

    parser.add_argument("SAM",
            nargs="?", type=argparse.FileType("r"), default=sys.stdin)

    return parser.parse_known_args()

def main():
    
    args, extra = parse_args()

    valid_spliced_reads = 0
    problem_reads = 0

    for ln in fileinput.input(extra):

        # Skip SAM header lines.
        if ln[0] == "@":
            continue

        ln_split = ln.split()

        # Get the CIGAR string.
        cigar = ln_split[5]

        # Look for N (skipped bases on the reference)
        if "N" in cigar:

            # Skip alignments with a low quality score.
            quality_score = int(ln_split[4])
            if quality_score < args.min_quality:
                continue
                    
            try:
                intron_length = int(cigar.split("N")[0].split("M")[-1])
                cs_split_M = cigar.split("M")
                edge5 = int(cs_split_M[0])
                edge3 = int(cs_split_M[-2].split("N")[-1])
            except ValueError:
                problem_reads += 1
                continue

            # Intron length must be > 50 and 6nt must map into each exon.
            # if intron_length > 50 and min_edge >= 6: 
            if intron_length >= args.min_intron_length \
                    and edge5 >= args.min_overlap\
                    and edge3 >= args.min_overlap:

                valid_spliced_reads += 1

                # Print a status as the script runs.
                if valid_spliced_reads % 100000 == 0:
                    sys.stderr.write(
                            "{} valid, {} problematic spliced reads\n"
                            .format(valid_spliced_reads, problem_reads))

                sys.stdout.write(ln)

if __name__ == "__main__":
    main()
