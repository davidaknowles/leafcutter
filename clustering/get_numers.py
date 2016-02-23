
import sys
import gzip


input_file=gzip.open(sys.argv[1], 'rb')
first_line=True
for l in input_file:
    if first_line:
        print " ".join(l.strip().split(" ")[1:]) # print the sample names
        first_line=False
    else:
        l=l.strip()
        words=l.split(" ")
        print(words[0] + " " + " ".join( [ g.split("/")[0] for g in words[1:] ] )  )

input_file.close()
