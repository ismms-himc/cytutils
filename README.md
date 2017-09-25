# cytutils
Cytometry quality control and reproducibility utilities. The methods available
in this package were developed by the Human Immune Monitoring Center at the
Icahn School of Medicine at Mount Sinai.

## Installation

TODO

## Available Methods

### Jensen-Shannon Divergence

The [Jensen-Shannon (JS) divergence](https://en.wikipedia.org/wiki/Jensen%E2%80%93Shannon_divergence)
is a statistical method for measuring the similarity between two probability
distributions. In the context of mass and flow cytometry, the JS divergence has
been utilized the compare the two-dimensional maps resulting from dimensionality
reduction methods such as PCA and t-SNE. Briefly, the maps are converted into
a probability distribution using kernel density estimation (KDE), which are in
then compared. Additional information is available in the supplementary
material of:

* Amir ED, Davis KL, Tadmor MD, Simonds EF, Levine JH, Bendall SC, Shenfeld DK,
Krishnaswamy S, Nolan GP, Pe'er D. viSNE enables visualization of high
dimensional single-cell data and reveals phenotypic heterogeneity of leukemia.
Nat Biotechnol. 2013 Jun;31(6):545-52.

`calculateJsDivergence` calculates the JS divergence between two probability
distributions. `calculateDrJsDivergence` calculates it between two matrices
which are the result of dimensionality reduction. For example:

```r
x <- read.csv("visne_map_1.csv")
y <- read.csv("visne_map_2.csv")

calculateDrJsDivergence(x, y)
```
