# leafcutter
LeafCutter: Annotation-free quantification of RNA splicing. 

Yang I Li, David A Knowles, Jonathan K Pritchard. 

For details please see
http://biorxiv.org/content/early/2016/03/16/044107

## Installation

To compile the R package to perform differential splicing analysis and make junction plots you can either...

### Option 1. Install from source
```
cd leafcutter
R CMD INSTALL --build .
```
You'll need the following R packages: `Rcpp, rstan, foreach, ggplot2, R.utils, gridExtra, reshape2, Hmisc, dplyr, doMC, optparse`. 

### Option 2. Install using devtools

This has the advantage of installing the required package dependencies for you. 
```
library(devtools)
install_github("davidaknowles/leafcutter/leafcutter")
```

## Usage

For a (hopefully) complete example of the complete pipeline, take a look at
```
example_data/worked_out_example.sh
```
:warning: the test data download is 4Gb! 

LeafCutter has two main components: 

1. Python code to 
   - generate intron excision counts from `junc` files (which can be obtained easily from `.bam` files)
   - group introns into clusters
2. `R` code to 
   * perform differential splicing (here, differential intron excision) analysis
   * plot (differentially spliced) clusters
   
:bug: we don't currently have scripts for splicing QTL calling in the git repo, please let us know if this would be useful for you! 

### Step 1. Converting `bam`s to `junc`s

I'm skipping Step 0 which would be mapping `fastq` files, e.g. using STAR, to obtain `.bam` files. 

We provide a helper script `scripts/bam2junc.sh` to (you guessed it) convert `bam` files to `junc` files. This step uses the CIGAR strings in the `bam` to quantify the usage of each intron. 

`example_data/worked_out_example.sh` gives you an example of how to do this in batch:
```
for bamfile in `ls my_awesome_data/*.bam`
do
    echo Converting $bamfile to $bamfile.junc
    sh ../scripts/bam2junc.sh $bamfile $bamfile.junc
    echo $bamfile.junc >> awesome_juncfiles.txt
done
```

This step is pretty fast but if you have samples numbering in the 100s you might want to do this on a cluster. Note that we also make a list of the generated `junc` files in `awesome_juncfiles.txt`. 

### Step 2. Intron clustering

Next we need to define intron clusters using the `leafcutter_cluster.py` script. For example: 

```
python ../clustering/leafcutter_cluster.py -j awesome_juncfiles.txt -m 50 -o awesome -l 500000
```

This will cluster together the introns fond in the `junc` files listed in `awesome_juncfiles.txt`, requiring 50 split reads supporting each cluster and allowing introns of up to 500kb. The predix `awesome` means the output will be called `awesome_perind_numers.counts.gz` (perind meaning these are the *per individual* counts). 

You can quickly check what's in that file with 
```
zcat awesome_perind_numers.counts.gz | more 
```
which should look something like this: 
```
TODO
```

### Step 3. Differential intron excision analysis

We can now use our nice intron count file to do differential splicing (DS) analysis, but first we need to make a file to specify which samples go in each group. In `worked_out_example.sh` there's some code to generate this file, but you can just make it yourself: it's just a two column tab-separated file where the first column is the sample name (i.e. the filename of the 'bam' without the extension) and the second is the column (what you call these is arbitrary, but note :bug: the command line interface currently only supports two groups). 

For the worked example this file, named `test_diff_introns.txt` looks like this: 
```
TODO
```

Having made that file we can run DS (this assumes you have successfully installed the `leafcutter` R package as described under Installation above) 
```
../scripts/leafcutter_ds.R --num_threads 4 ../example_data/testYRIvsEU_perind_numers.counts.gz ../example_data/test_diff_intron.txt
```

Running `../scripts/leafcutter_ds.R -h` will give usage info for this script. 

Two tab-separated text files are output:

1. `leafcutter_ds_cluster_significance.txt`. This shows per cluster `p`-values for there being differential intron excision between the two groups tested. The columns are
 1. cluster: the cluster id
 2. Status: whether this cluster was a) successfully tested b) not tested for some reason (e.g. too many introns) c) there was an error during testing - this should be rare. 
 3. loglr: log likelihood ratio between the null model (no difference between the groups) and alternative (there is a difference) 
 4. df: degrees of freedom, equal to the number of introns in the cluster minus one (assuming two groups)
 5. p: the resulting p-value under the asymptotic Chi-squared distribution

2. `leafcutter_ds_effect_sizes.txt`. This shows per intron effect sizes between the groups, with columns:
 1. intron: this has the form chromosome:intron_start:intron_end:cluster_id
 2. effect size. 

### Step 4. Plotting splice junctions

This will make a pdf with plots of the differentially spliced clusters detected at an FDR of 5%. 
```
../scripts/ds_plots.R -e ../leafcutter/data/gencode19_exons.txt.gz ../example_data/awesome_perind_numers.counts.gz ../example_data/test_diff_intron.txt leafcutter_ds_cluster_significance.txt -f 0.05
```