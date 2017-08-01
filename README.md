# leafcutter

<img src="./logo.png" width="200"> **Annotation-free quantification of RNA splicing.**

*Yang I Li, David A Knowles, Jonathan K Pritchard.*

Leafcutter quantifies RNA splicing variation using short-read RNA-seq data. The core idea is to leverage spliced reads (reads that span an intron) to quantify (differential) intron usage across samples. The advantages of this approach include
* easy detection of novel introns
* modeling of more complex splicing events than exonic PSI
* avoiding the challenge of isoform abundance estimation
* simple, computationally efficient algorithms scaling to 100s or even 1000s of samples

For details please see our [bioRxiv preprint](http://biorxiv.org/content/early/2016/03/16/044107)

Full documentation is available at <http://davidaknowles.github.io/leafcutter/>
