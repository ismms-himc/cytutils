# cytutils
Cytometry quality control and reproducibility utilities. The methods available
in this package were developed by the Human Immune Monitoring Center at the
Icahn School of Medicine at Mount Sinai. Please contact Adeeb Rahman at
adeeb.rahman@mssm.edu with suggestions and questions.

## Installation

The `cytutils` package is currently only available as a direct installation from
github. You can use the `devtools` package to install it:

```r
install.packages("devtools")
library(devtools)
install_github("ismmshimc/cytutils")
```

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
[PubMedLink](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4076922/)

`generate2dJsDivergenceDataFrame` accepts a list of source file paths as an argument 
and calculates the JS divergence between every possible pair combination of
two-dimensional matrices. The two-dimensional matrices are generated from the source 
files and columns of interest (`two_d_map_column_names`). The matrices may or 
may not be a result of dimensionality reduction. By default, 
`two_d_map_column_names` is set to Cytobank's default tSNE column names.

```r
source_filepaths <- c("/path/to/abc.fcs", "/path/to/def.fcs", "/path/to/ghi.fcs")
generate2dJsDivergenceDataFrame(source_filepaths, two_d_map_column_names = c("tSNE1", "tSNE2"))
```

The function returns a data frame like the below:

| file1 | file2 | divergence_value |
| ------ | ------ |------ |
| abc.fcs | def.fcs | 0.010899564 |
| abc.fcs | ghi.fcs | 0.005387104 |
| def.fcs | ghi.fcs | 0.051054440 |


### Average Overlap Frequency

The Average Overlap Frequency (AOF) is a simple, intuitive metric for the
staining quality of cytometry data. The method assumes that channels have a
bimodal distribution, corresponding to negative and positive populations. It
then quantifies how many negative cells are "positive-looking" and vice versa.
The AOF returns a value between 0 and 1, where 0 signifies perfect separation.
An in-depth explanation can be found here:

* Amir ED, Guo XV, Mayovska O, Rahman A. Average Overlap Frequency: A simple
metric to evaluate staining quality and community identification in high
dimensional mass cytometry experiments. J Immunol Methods. 2017 Sep 4.
[PubMed Link](https://www.ncbi.nlm.nih.gov/pubmed/28882613)

If the negative and positive populations for a given marker are known (for
example, using manual gating) use the `calculateAof` function. Alternatively,
we include the greedy search algorithm that uses a clustering of the data to
find (for each marker) a partition that minimizes the AOF. The greedy method is
available using the `greedyCytometryAof` function. For example:

```r
fcs_data <- flowCore::read.FCS("fcs_data.fcs")
# Load output of semi-supervised clustering method.
y <- read.csv("clustering.csv")

# Calculate AOF for CD3 using T cells versus non-T cells.
x <- fcs_data@exprs[, "CD3"]
t_cell_indices <- grep("t_cell", y)
calculateAof(x, t_cell_indices)

# Use greedy search to optimize AOF for each channel.
channel_names <- c("CD3", "CD19", "CD33", "CD56", "CD42")
greedyCytometryAof(fcs_data@exprs, y, channel_names)
```
