

def main(flist, outPrefix):
    
    N = 80

    lsts = []
    for ln in open(flist):
        lsts.append(ln.strip())
    print len(lsts)
    tmpfiles = []
    while len(lsts) > 1:    
        clst = []
        
        for i in range(0,(len(lsts)/N)+1): 
            lst = lsts[N*i:N*(i+1)]
            if len(lst) > 0:
                clst.append(lst)
        lsts = []
        print [len(x) for x in clst]
        for lst in clst:
            if len(lst) == 0: continue
            tmpfile = tempfile.mktemp()
            os.mkdir(tmpfile)
            foutname = tmpfile+"/tmpmerge.gz"
            fout = gzip.open(foutname,'w')
            print foutname
            merge_files(lst, fout)
            lsts.append(foutname)
            tmpfiles.append(foutname)
            fout.close()
    
    shutil.move(lsts[0], outPrefix+"_perind.counts.gz")
    #for tmpfile in tmpfiles:
    #    try: os.remove(tmpfile)
    #    except: continue

def merge_files(fnames, fout):
    fopen = []
    for fname in fnames:
        if fname[-3:] == ".gz":
            fopen.append(gzip.open(fname))
        else:
            fopen.append(open(fname))

    finished = False
    N = 0
    while not finished:
        N += 1
        if N % 50000 == 0: 
            sys.stderr.write("read %d...\n"%N)
        buf = []
        for f in fopen:
            ln = f.readline().split()
            if len(ln) == 0: 
                finished = True
                break
            chrom = ln[0]
            data = ln[1:]
            if len(buf) == 0:
                buf.append(chrom)
            buf += data
        if len(buf) > 0:
            if buf[0] == "chrom":
                print buf[:3], len(buf)
            fout.write(" ".join(buf)+'\n')
        else:
            break

    
    for fin in fopen:
        fin.close()

if __name__ == "__main__":

    import sys
    import tempfile
    import os
    import gzip
    import shutil
    outPrefix = sys.argv[1]
    filelist = outPrefix+"_sortedlibs"
    main(filelist, outPrefix)
