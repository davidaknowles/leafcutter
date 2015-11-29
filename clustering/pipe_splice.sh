juncFiles="../run/intron_files/*.junc"
outPrefix="../run/GEUVADIS"

maxIntronLen=100000
minReadsClu=50

echo "Checking junction files..."
python check_files.py $outPrefix\_chrom $juncFiles

echo "Pooling individuals to identify intron clusters..."
python get_intron_splicing_pooled.py $outPrefix $maxIntronLen $juncFiles

echo "Spliting intron clusters..."
python refine_intron_splicing.py $outPrefix $minReadsClu

echo "Sorting junction files..."
python sort_junc.py $outPrefix $juncFiles   

echo "Getting per individuals counts..."
python merge_junc.py $outPrefix

echo "Assigning genes to clusters..."
python get_cluster_gene.py $outPrefix\_perind.counts.gz
