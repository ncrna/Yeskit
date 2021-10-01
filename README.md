# Yeskit (*Yet another single-cell analysis toolkit*) <img src="https://github.com/ncrna/Yeskit/blob/master/inst/figures/Yeskit_logo.png" align="right" height=142 width=164/>

`Yeskit` is an R package designed for single-cell gene expression data _importation_, _integration_, _clustering_, _differential analysis_, _functional analysis_, and _visualization_. Since `Yeskit` does not change the default data structure of `Seurat`, it can be easily integrated into most existing scRNA-seq analysis workflows. 

`Yeskit` can be used to read other information (such as gene mutation-by-cell matrix, pathogen count-by-cell matrix) and store them as additional data in the Seurat obj@metadata slot. In addition, `Yeskit` also has the function of reading and distinguishing the source of cells from the scRNA-seq data of xenografts samples (PDX model).

When there are many points in the vector diagram, editing becomes difficult. To this end, most visualization functions in `Yeskit` have the option to rasterize the _geom_point_ layer of ggplot2 and keep all axes, labels, and text in vector format.

## Installation

`Yeskit` package can be easily installed under the following instruction:

```
devtools::install_github("ncrna/Yeskit")
```

## Vignettes

Please check the following link to get a rapid understanding about `Yeskit`:

[Walkthrough - Decoding Intracellular Pathogens of H3N2 at the Single-Cell level using Yeskit](https://htmlpreview.github.io/?https://github.com/ncrna/Yeskit/blob/main/vignettes/Decoding_Intracellular_Pathogens_of_H3N2_at_the_Single-Cell_level_using_Yeskit.html)

## Troubleshooting

In a particular case, the Cairo package can be loaded but crashes the R process when called. This is probably the case that the version of Cairo was built with an older version of R. Try the following instruction:
```
Bioconductor::install("Cairo", force = TRUE)
```

For OS X users, you may happen to an error like this:
```
Error : .onLoad failed in loadNamespace() for ‘Cairo’
... 
Library not loaded: /opt/X11/lib/libXrender.1.dylib
```
The reason is X11 doesn't ship with OS X any more, users can download and install **XQuartz** from https://www.xquartz.org to solve this problem.

## Help

If you have any problems, comments or suggestions, please feel free to contact _**Wei Zhang**_, [admin@ncrna.net](mailto:admin@ncrna.net).

## How to cite

**Decoding Intracellular Pathogens of scRNA-seq experiments with PathogenTrack and Yeskit.**








