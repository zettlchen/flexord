# flexord

This `R` package and `github` repository is an add-on package to the
`R` packages [`flexclust`](https://cran.r-project.org/web/packages/flexclust/index.html)
and [`flexmix`](https://cran.r-project.org/web/packages/flexmix/index.html). These two root
packages are suites for flexible method comparison and clustering with both partitioning
and model-based methods.

With `flexord`, we provide new distance and centroid functions and new distribution drivers
that are tailored towards clustering ordinal or mixed-with-ordinal data.


We wrote this package accompanying the paper "Ordinal clustering in the flex-scheme" (2025)
by Ernst D, Ortega Menjivar L, Scharl T, Grün B, currently submitted to [AJS](https://ajs.or.at).
In this paper, we collected the state of the art of ordinal data clustering and applied many
of the methods in our extensive simulation study. The simulation reproduction code itself is collected
in the repository [`AJS-flexord`](https://github.com/zettlchen/AJS-flexord). <replace with Zenodo DOI>


You can find a more end user friendly introduction into the package, and a listing of the new
methods we provide [here](https://dernst.github.io/flexord/articles/Intro2Flexord.html).

## Installation

You can install the developer version of the package via:

```
devtools::install_github("dernst/flexord")
```

The stable version will soon be available on CRAN (via `install.packages("flexord")`).
