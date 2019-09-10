#!/bin/bash

# Manually working through the differential splicing pipeline.

# wget bam files (4Gb!)
wget https://www.dropbox.com/s/pni1zq5y6cr4tx5/example_geuvadis.tar.gz?dl=0 -O example_geuvadis.tar.gz
tar -xvf example_geuvadis.tar.gz

# Convert bam to junction files
if [ -e test_juncfiles.txt ]; then rm test_juncfiles.txt; fi

for bamfile in `ls example_geuvadis/*chr1.bam`; do
    echo Converting $bamfile to $bamfile.junc
    sh ../scripts/bam2junc.sh $bamfile $bamfile.junc
    echo $bamfile.junc >> test_juncfiles.txt
done

# Finds intron clusters and quantifies junction usage within them.
python ../clustering/leafcutter_cluster.py -j test_juncfiles.txt -m 50 -o testYRIvsEU -l 500000

# Differential splicing analysis.
../scripts/leafcutter_ds.R --num_threads 4 ../example_data/testYRIvsEU_perind_numers.counts.gz example_geuvadis/groups_file.txt

# Plot pdf differentially spliced clusters.
../scripts/ds_plots.R -e ../leafcutter/data/gencode19_exons.txt.gz ../example_data/testYRIvsEU_perind_numers.counts.gz example_geuvadis/groups_file.txt leafcutter_ds_cluster_significance.txt -f 0.05

## Leafviz visualization shiny app
pushd ../leafviz/
# Download hg19 annotation files
./download_human_annotation_codes.sh
# Cache results for viewing
./prepare_results.R --meta_data_file ../example_data/example_geuvadis/groups_file.txt \
  --code leafcutter ../example_data/testYRIvsEU_perind_numers.counts.gz \
  ../example_data/leafcutter_ds_cluster_significance.txt \
  ../example_data/leafcutter_ds_effect_sizes.txt \
  annotation_codes/gencode_hg19/gencode_hg19 \
  -o testYRIvsEU.RData
# Launch shiny app locally
./run_leafviz.R testYRIvsEU.RData
popd
