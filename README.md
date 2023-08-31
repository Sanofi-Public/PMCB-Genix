[![Project Status](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![Lifecycle](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://www.tidyverse.org/lifecycle/#stable)

<p align="center" width="100%">
<img width="30%" src="vignettes/genix-logo.png"> 
</p>

# Sargent 

Sargent is a transformation- and cluster-free cell-type annotation method that 
operates at individual cell resolution by applying a scoring system to scRNA-seq 
data based on sets of marker genes associated with cell types.


## Installation

**sargent** can be installed directly from this github with:

```{r}
if (!require("devtools", quietly = TRUE))
  install.packages("devtools")

devtools::install_github("Sanofi-Public/PMCB-Sargent", 
                         build_vignettes=FALSE)
```

If you wish to build a local version of the vignette use:

```{r}
if (!require("Seurat", quietly = TRUE))
  install.packages("Seurat")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("BiocStyle")

devtools::install_github("Sanofi-Public/PMCB-Sargent", 
                         build_vignettes=TRUE)
```


## Getting started

Once installed the best place to get started is the [vignette][vignette].


## Contact

For help and questions please contact the [sargent's maintenance team](mailto:nima.nouri@sanofi.com).


## Citing Sargent

If you use Sargent please cite our paper: Nouri N. et al. "A marker gene-based 
method for identifying the cell-type of origin from single-cell RNA sequencing 
data", MethodsX, 2023, [doi:10.1016/j.mex.2023.102196][paper].

```
  @Article{,
    author = {Nima Nouri et al.},
    title = {A marker gene-based method for identifying the cell-type of origin 
             from single-cell RNA sequencing data},
    journal = {MethodsX},
    year = {2023},
    url = {https://www.sciencedirect.com/science/article/pii/S2215016123001966},
    doi = {10.1016/j.mex.2023.102196},
  }
```

[vignette]: https://github.com/Sanofi-Public/PMCB-Sargent/blob/master/vignettes/Sargent-Vignette.Rmd
[bioc]: https://bioconductor.org/packages/devel/bioc/html/sargent.html
[paper]: https://www.sciencedirect.com/science/article/pii/S2215016123001966
