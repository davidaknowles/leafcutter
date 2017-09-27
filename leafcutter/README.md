# LeafCutter: Annotation-free quantification of RNA splicing

<img src="./logo.png" width="200"> 

[Yang I. Li](http://web.stanford.edu/~yangili/index.html)<sup>1</sup>, [David A. Knowles](http://cs.stanford.edu/people/davidknowles/)<sup>1</sup>, Jack Humphrey, Alvaro N. Barbeira, Scott P. Dickinson, Hae Kyung Im, [Jonathan K. Pritchard](http://web.stanford.edu/group/pritchardlab/home.html)

<sup>1</sup>*Equal contribution*

Leafcutter quantifies RNA splicing variation using short-read RNA-seq data. The core idea is to leverage spliced reads (reads that span an intron) to quantify (differential) intron usage across samples. The advantages of this approach include
* easy detection of novel introns
* modeling of more complex splicing events than exonic PSI
* avoiding the challenge of isoform abundance estimation
* simple, computationally efficient algorithms scaling to 100s or even 1000s of samples

For details please see our [bioRxiv preprint](http://www.biorxiv.org/content/early/2017/09/07/044107)

* [Installation](./articles/Installation.html)
* [Differential splicing](./articles/Installation.html)
* [Visualization](./articles/Visualization.html)
* [SplicingQTL](./articles/sQTL.html)

Check out a demo leafcutter [shiny](https://shiny.rstudio.com/) app [here](https://leafcutter.shinyapps.io/leafviz/): 10 brain vs. 10 heart samples from [GTEx](https://www.gtexportal.org/home/). 

We have a Google group for user questions at <https://groups.google.com/forum/#!forum/leafcutter-users>
