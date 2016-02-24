
leafCutterDir='.'
bamfile=$1

samtools view $bamfile | python $leafcCutterDir/scripts/filter_cs.py | $leafcCutterDir/scripts/sam2bed.pl --use-RNA-strand - $s.junc