# LeafCutter: Annotation-free quantification of RNA splicing

<img src="./logo.png" width="200"> 

[Yang I. Li](https://ggsb.uchicago.edu/program/faculty/yang-li)<sup>1</sup>, [David A. Knowles](https://daklab.github.io/)<sup>1</sup>, Jack Humphrey, Alvaro N. Barbeira, Scott P. Dickinson, Hae Kyung Im, [Jonathan K. Pritchard](http://web.stanford.edu/group/pritchardlab/home.html)

<sup>1</sup>*Equal contribution*

Leafcutter quantifies RNA splicing variation using short-read RNA-seq data. The core idea is to leverage spliced reads (reads that span an intron) to quantify (differential) intron usage across samples. The advantages of this approach include
* easy detection of novel introns
* modeling of more complex splicing events than exonic PSI
* avoiding the challenge of isoform abundance estimation
* simple, computationally efficient algorithms scaling to 100s or even 1000s of samples

For details please see our [bioRxiv preprint](http://www.biorxiv.org/content/early/2017/09/07/044107) and corresponding [Nature Genetics publication](https://www.nature.com/articles/s41588-017-0004-9).

Additionally, for full details on the leafcutter for Mendelian Diseases (leafcutterMD) method that performs outlier splicing detection, see our [Bioinformatics publication](http://dx.doi.org/10.1093/bioinformatics/btaa259).

* [Installation](./articles/Installation.html)
* [Differential splicing](./articles/Usage.html)
* [Outlier splicing](./articles/UsageLeafcutterMD.html)
* [Visualization](./articles/Visualization.html)
* [SplicingQTL](./articles/sQTL.html)

Check out a demo leafcutter [shiny](https://shiny.rstudio.com/) app [here](https://leafcutter.shinyapps.io/leafviz/): 10 brain vs. 10 heart samples from [GTEx](https://www.gtexportal.org/home/). 

We have a Google group for user questions at <https://groups.google.com/forum/#!forum/leafcutter-users>
