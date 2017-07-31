
leafCutterDir='/home/yangili/tools/leafcutter/'
bamfile=$1
bedfile=$1.bed
juncfile=$2
samtools view $bamfile | python $leafCutterDir/scripts/filter_cs.py | $leafCutterDir/scripts/sam2bed.pl --use-RNA-strand - $bedfile
$leafCutterDir/scripts/bed2junc.pl $bedfile $juncfile
rm $bedfile
