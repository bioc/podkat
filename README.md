# PODKAT - An R Package for Association Testing Involving Rare and Private Variants
This package provides an association test that is capable of dealing with very rare and even private variants. This is accomplished by a kernel-based approach that takes the positions of the variants into account. The test can be used for pre-processed matrix data, but also directly for variant data stored in VCF files. Association testing can be performed whole-genome, whole-exome, or restricted to pre-defined regions of interest. The test is complemented by tools for analyzing and visualizing the results.

## Installation

The package can be installed from
[Bioconductor](https://bioconductor.org/). Therefore, the the simplest way to install the package is to enter
```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("podkat")
```
into your R session. If, for what reason ever, you prefer to install the package manually, follow the instructions in the [user manual](https://bioconductor.org/packages/release/bioc/vignettes/podkat/inst/doc/podkat.pdf).

## User support

If you encounter any issues or if you have any question that might be of interest also for other users, before writing a private message to the package developers/maintainers, please create an issue in this repository and also consider posting on [Bioconductor Support](https://support.bioconductor.org/) or on [StackOverflow](https://stackoverflow.com/). For other matters regarding the package, please contact the package author.

## Citing this package

If you use this package for research that is published later, you are kindly asked to cite it as follows:

- U. Bodenhofer. PODKAT: An R Package for Association Testing Involving Rare and Private Variants. R package.
  DOI: [10.18129/B9.bioc.podkat](https://doi.org/10.18129/B9.bioc.podkat)
