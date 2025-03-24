# flexord

This `R` package and `github` repository is an add-on package to the
`R` packages [`flexclust`](https://CRAN.R-project.org/package=flexclust/)
and [`flexmix`](https://CRAN.R-project.org/package=flexmix/). These two root
packages are suites for flexible clustering with both partitioning
and model-based methods allowing for easy method variation and comparison.

With `flexord`, we provide new distance and centroid functions and new
model drivers for component distributions that are tailored towards
clustering ordinal or mixed-with-ordinal data.

We wrote this package to accompany the paper "Ordinal clustering in the flex-scheme" (2025)
by Ernst D, Ortega Menjivar L, Scharl T, Grün B, currently submitted to [AJS](https://ajs.or.at).
In this paper, we reviewed methods available for ordinal data clustering and applied many
of the methods in our extensive simulation study. The replication code for the simulation study itself is collected
in the repository [`AJS-flexord`](https://doi.org/10.5281/zenodo.15074618).

You can find a more end user friendly introduction to the package, and a listing of the new
methods we provide [here](https://dernst.github.io/flexord/articles/Intro2Flexord.html).

## Installation

The stable version of the package is available at CRAN, and can be installed via:
```
install.packages("flexord")
```

You can install the developer version of the package via:

```
devtools::install_github("dernst/flexord")
```
