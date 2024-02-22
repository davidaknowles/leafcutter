# leafcutter

[Yang I. Li](https://thelilab.com/)<sup>1</sup>, [David A. Knowles](https://daklab.github.io/)<sup>1</sup>, [Jack Humphrey](https://jackhump.github.io/), Alvaro N. Barbeira, Scott P. Dickinson, Hae Kyung Im, [Jonathan K. Pritchard](http://web.stanford.edu/group/pritchardlab/home.html)

This branch has a Python/pytorch/pyro reimplementation of the LeafCutter differential splicing test. Users (including us!) often had problems installing the RStan based package. 

You should be able to install with `python -m pip install leafcutter`. 

The installation process will install `leafcutter-cluster` and `leafcutter-ds` as command line tools. 

<img src="https://github.com/davidaknowles/leafcutter/blob/master/docs/logo.png" width="200"> **Annotation-free quantification of RNA splicing.**

Leafcutter quantifies RNA splicing variation using short-read RNA-seq data. The core idea is to leverage spliced reads (reads that span an intron) to quantify (differential) intron usage across samples. The advantages of this approach include
* easy detection of novel introns
* modeling of more complex splicing events than exonic PSI
* avoiding the challenge of isoform abundance estimation
* simple, computationally efficient algorithms scaling to 100s or even 1000s of samples

For details please see our [bioRxiv preprint](http://www.biorxiv.org/content/early/2017/09/07/044107) and corresponding [Nature Genetics publication](https://www.nature.com/articles/s41588-017-0004-9).

Additionally, for full details on the leafcutter for Mendelian Diseases (leafcutterMD) method that performs outlier splicing detection, see our [Bioinformatics publication](http://dx.doi.org/10.1093/bioinformatics/btaa259).

Full documentation is available at <http://davidaknowles.github.io/leafcutter/>

If you have usage questions we've setup a Google group here: <https://groups.google.com/forum/#!forum/leafcutter-users>

We've developed a leafcutter [shiny](https://shiny.rstudio.com/) app for visualizing leafcutter results: you can view an example [here](https://leafcutter.shinyapps.io/leafviz/). This shows leafcutter differential splicing results for a comparison of 10 brain vs. 10 heart samples (5 male, 5 female in each group) from [GTEx](https://www.gtexportal.org/home/). 
