
juncFilesDir="~/scailscratch/splicing/processed/*.junc"
juncFiles="junc_files.txt"
outPrefix="~/scailscratch/splicing/gtex_final"

#juncFiles="../run/intron_files/*.junc"
#outPrefix="../run/GEUVADIS"

mkdir outPrefix

maxIntronLen=100000
minReadsClu=50

echo "Checking junction files..."
#eval "ls $juncFilesDir > $juncFiles" 
#eval "python check_files.py $outPrefix\_chrom $juncFiles"

echo "Sorting junction files..."
#eval "python sort_junc.py $outPrefix $juncFiles"

echo "Pooling individuals to identify intron clusters..."
#eval "python get_intron_splicing_pooled.py $outPrefix $maxIntronLen $juncFiles"

echo "Spliting intron clusters..."
#eval "python refine_intron_splicing.py $outPrefix $minReadsClu"

echo "Getting per individuals counts..."
#eval "python get_intron_splicing_perind.py $outPrefix"

#eval "python merge_junc.py $outPrefix"

echo "Assigning genes to clusters..."
eval "python get_cluster_gene.py $outPrefix\_perind.counts.gz"
