#!/usr/bin/env python

import sys

valid_spliced_reads=0
problem_reads=0

for ln in sys.stdin:
    if ln[0] == "@": continue
    ln_split=ln.split()
    cs = ln_split[5]
    if "N" in cs:
        quality_score=int(ln_split[4])
        if quality_score < 10: continue # quality score
                
        try:
            L = int(cs.split("N")[0].split("M")[-1])
            cs_split_M=cs.split("M")
            edge5, edge3 = cs_split_M[0], cs_split_M[-2].split("N")[-1]
            minedge = min([int(edge5), int(edge3)])
        except ValueError:
            problem_reads += 1
            continue
            
        if L > 50 and minedge >= 6: # intron length must be > 50 and 6nt must map into each exon
            valid_spliced_reads += 1
            if valid_spliced_reads % 100000 == 0: sys.stderr.write("%d valid, %d problematic spliced reads\n" % (valid_spliced_reads, problem_reads) )
            sys.stdout.write(ln)
        elif minedge < 6:
            pass
