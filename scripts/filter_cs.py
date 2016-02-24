import sys

N = 0
for ln in sys.stdin:
    if ln[0] == "@": continue
    cs = ln.split()[5]
    if "N" not in cs:
        pass
    else:
        if int(ln.split()[4]) < 10: continue # quality score
        L = int(cs.split("N")[0].split("M")[-1])
        
        edge5, edge3 = cs.split("M")[0], cs.split("M")[-2].split("N")[-1]
        try: minedge = min([int(edge5), int(edge3)])
        except: continue
        if L > 50 and minedge >= 6: # intron length must be > 50 and 6nt must map into each exon
            sys.stdout.write(ln)
        elif minedge < 6:
            pass
