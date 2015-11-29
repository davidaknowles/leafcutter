

def main(flist, outPrefix):
    

    lsts = []
    for ln in open(flist):
        lsts.append(ln.strip())
    sys.stderr.write("merging %d junction files...\n"%(len(lsts)))
    
    # Change 500 if max open file is < 500
    N = min([500, max([100, int(len(lsts)**(0.5))])])

    tmpfiles = []
    while len(lsts) > 1:    
        clst = []
        
        for i in range(0,(len(lsts)/N)+1): 
            lst = lsts[N*i:N*(i+1)]
            if len(lst) > 0:
                clst.append(lst)
        lsts = []
    
        for lst in clst:
            if len(lst) == 0: continue
            tmpfile = tempfile.mktemp()
            os.mkdir(tmpfile)
            foutname = tmpfile+"/tmpmerge.gz"
            fout = gzip.open(foutname,'w')
            
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
            sys.stderr.write(".")
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
                sys.stderr.write("merging %d files"%(len(buf)-1))
            fout.write(" ".join(buf)+'\n')
        else:
            break

    sys.stderr.write("done.\n")
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
