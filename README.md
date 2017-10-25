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
files and columns of interest (`two_d_map_column_names`). While the matrices may be 
the result of dimensionality reduction this is not required. By default, 
`two_d_map_column_names` is set to Cytobank's default tSNE column names.

```r
source_filepaths <- c("/path/to/A.fcs", "/path/to/B.fcs", "/path/to/C.fcs")
generate2dJsDivergenceDataFrame(source_filepaths, two_d_map_column_names = c("tSNE1", "tSNE2"))
```

The function returns a data frame like the below:

| file1 | file2 | divergence_value |
| ------ | ------ |------ |
| A.fcs | B.fcs | 0.010899564 |
| A.fcs | C.fcs | 0.005387104 |
| B.fcs | C.fcs | 0.051054440 |


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


#### Generalized workflow:
```r
fcs_data <- flowCore::read.FCS("fcs_data.fcs")
# Load output of semi-supervised clustering method. y is a vector of assignments.
y <- read.csv("clustering.csv")

# Use greedy search to optimize AOF for each channel.

channel_names <- c("CD3", "CD19", "CD33", "CD56", "CD42")
greedyCytometryAof(fcs_data@exprs, y, channel_names)

# Calculate AOF for Er168Di (CD3) using T cells versus non-T cells.
x <- fcs_data@exprs[, "Er168Di"]
t_cell_indices <- grep("t_cell", y)
non_t_cell_indices <- grep("t_cell", y, invert = TRUE)
# Note: Entering a cofactor will result in data being transformed. If x has not 
# yet been transformed, a cofactor should be provided to the function below.
calculateAof(x, t_cell_indices, non_t_cell_indices, cofactor = 5)
```


#### A possible workflow for when using manually gated data:
```r
# Use manually gated data to assign sample cells to specific populations
sample_1_base_fcs_data <- flowCore::read.FCS("sample_1_base.fcs")
manual_labeling_filepath <- "/path/to/samples_manual_labeling.csv"
samples_filepath <- "/path/to/samples.csv"
data_dir <- "/path/to/fcs_files"

single_sample_labels <- generatePopulationAssignments(
							manual_labeling_filepath, 
							samples_filepath, 
							data_dir)

# Use greedy search to optimize AOF for each channel.
channel_names <- c("Er168Di", "Nd142Di", "Gd158Di", "Dy161Di")

# Load output of semi-supervised clustering method and generate a vector of assignments
# (cell_assignments_ordered below).
cell_assignments <- read.csv("cell_assignments.csv")
index_label_pairs <- cell_assignments[c("Index", "Label")]
expected_num_rows <- nrow(sample_1_base_fcs_data@exprs)

# We find the indices that were not labeled and add them to our "Label" column 
# with a label of NA
unlabeled_indices <- setdiff(1:expected_num_rows, index_label_pairs$Index)
na_assignments <- rep(NA, length(unlabeled_indices))
unlabeled_assignments <- data.frame(unlabeled_indices, na_assignments)
colnames(unlabeled_assignments) <- c("Index", "Label")

index_label_pairs_complete <- rbind(index_label_pairs, unlabeled_assignments)
index_label_pairs_ordered <- index_label_pairs_complete[order(index_label_pairs_complete$Index),]
cell_assignments_ordered <- as.vector(index_label_pairs_ordered$Label)

greedyCytometryAof(sample_1_base_fcs_data@exprs, cell_assignments_ordered, channel_names, cofactor = 5) # =>

#  ChannelName        Aof
# 1     Er168Di 0.050136740
# 2     Nd142Di 0.007091019
# 3     Gd158Di 0.32359932
# 4     Dy161Di 0.73802856

# Calculate AOF for Er168Di (CD3) using T cells versus non-T cells.
x <- sample_1_base_fcs_data@exprs[, "Er168Di"]

y <- single_sample_labels["sample_1"] # Note: "sample_1" was a sample_id in our samples csv file.
t_cell_indices <- grep(TRUE, y$sample_1$t_cell)
non_t_cell_indices <- grep(FALSE, y$sample_1$t_cell)
calculateAof(x, t_cell_indices, non_t_cell_indices, cofactor = 5) # =>  0.003321323
```

#### Calculating AOF for multiple channels.
```r
# Use manually gated data to assign sample cells to specific populations
sample_1_base_fcs_data <- flowCore::read.FCS("sample_1_base.fcs")
manual_labeling_filepath <- "/path/to/samples_manual_labeling.csv"
samples_filepath <- "/path/to/samples.csv"
data_dir <- "/path/to/fcs_files"

single_sample_labels <- generatePopulationAssignments(
							manual_labeling_filepath, 
							samples_filepath, 
							data_dir)

# Calculate AOF for Er168Di (CD3) and Nd142Di, designating positive and negative 
# populations for each channel.

sample_1_base_fcs_data <- flowCore::read.FCS("sample_1_base.fcs")
# The below is a csv with the first column representing channel names (i.e. Er168Di,
# Nd142Di, etc.). Column names start with "channel", followed by cell population
# names (i.e. "b_cell", "t_cell", etc.). Cell values are TRUE, FALSE, or left empty. 
channel_population_relationships_filepath <- "channel_population_relationships.csv"
sample_id <- "sample_1"
base_fcs_data_filepath <- "/path/to/base_fcs_data.fcs"

calculateMultiChannelAof(channel_population_relationships_filepath, 
						base_fcs_data_filepath, 
						single_sample_labels, 
						sample_id,
						cofactor = 5)

#   ChannelName       Aof
# 1     Er168Di 0.050136740
# 2     Nd142Di 0.007091019

```