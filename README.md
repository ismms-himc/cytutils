# cytutils
Cytometry quality control and reproducibility utilities. The methods available
in this package were developed by the Human Immune Monitoring Center at the
Icahn School of Medicine at Mount Sinai. Please contact Adeeb Rahman at
adeeb.rahman@mssm.edu with suggestions and questions.

## Installation

Install `libcurl-dev`. E.g., on Ubuntu type:

```bash
sudo apt-get install libcurl4-openssl-dev
```

The `cytutils` package is currently only available as a direct installation from
github. You can use the `devtools` package to install it:

```r
install.packages("devtools")
library(devtools)
install_github("ismmshimc/cytutils")
```

## Available Methods

### FCS Channel Renaming

The FCS channel renaming script allows you to easily rename channel names or
descriptions over a large set of FCS files. The main motivation are instances
where some files have mismatching descriptions. For example, naming a channel
CD8 in some samples and CD8a in others, or Nd142 (marker channel) and Ce142
(bead channel).

The current version of the method is a two-step process. We are going to use the
data from [Lavin et al., 2017](https://www.ncbi.nlm.nih.gov/pubmed/28475900) as
an example. If you would like to follow this example, you are welcome to
[download the data from FlowRepository](https://flowrepository.org/id/FR-FCM-ZY9N).
We assume that the FCS files have been unzipped to `D:/data/lavin`.

Once you have the FCS files unzipped, first run the script to generate a
`channel_rename.csv` CSV file:

```{r}
> cytutils::channelRename("D:/data/lavin")
Generating a new channel_rename.csv file
	D:/data/lavin/245_myeloid_blood.fcs
	D:/data/lavin/245_myeloid_lung.fcs
	D:/data/lavin/245_myeloid_tumor.fcs
	.
	.
	.
channel_rename.csv created. Update the file and re-run channelRename to export new FCS files
```

This will result in a [`channel_rename.csv`](examples/channel_rename_original.csv)
file, which includes channel information from all of the FCS files in the target
directory. In other words, each unique channel (combination of mass, name, and
description) will have a row in this file.

The file includes six columns. The first three (`mass`, `name`, and `desc`) are
the existing channels. The next column, `dup`, is `TRUE` if this mass is 
duplicated across files -- at least two files have this mass, but the respective
channels have different name or description. This column is supplied to help you
find channels that might have been named differently, and does not affect the
operation of the follow-up step.

The final two columns, `new_name` and `new_desc`, are initially identical to
`name` and `desc`. These are the columns that will be used for renaming. Set the
`new_name` and `new_desc` according to the new values you want. Download
[`channel_rename.csv`](examples/channel_rename.csv) for one suggestion. Then,
run the script again (do not rename `channel_rename.csv`!):

```{r}
> cytutils::channelRename("D:/data/lavin")
Importing channel_rename.csv and exporting new FCS files
	D:/data/lavin/245_myeloid_blood.fcs
	--> D:/data/lavin/channel_rename/245_myeloid_blood.fcs.renamed.fcs
	D:/data/lavin/245_myeloid_lung.fcs
	--> D:/data/lavin/channel_rename/245_myeloid_lung.fcs.renamed.fcs
	D:/data/lavin/245_myeloid_tumor.fcs
	--> D:/data/lavin/channel_rename/245_myeloid_tumor.fcs.renamed.fcs
	.
	.
	.
Export done
There were 50 or more warnings (use warnings() to see the first 50)
```

*(Please disregard the warnings. They are generated due to flowCore's
`write.FCS` still being considered "experimental".)*

The script will import each FCS file, rename the relevant channels, and export
a new FCS file under the `/channel_rename/` directory. In order to avoid
accidentally overwriting the original FCS files, the new files will have the
`.renamed.fcs` suffix added to them.

We would like to thank cytoforum for the [thread that inspired this feature](http://cytoforum.stanford.edu/viewtopic.php?f=3&t=874#p2536).

> The package also includes `importChannelNames` and `renameFcsFileChannels`
> which decouple the two steps and provide the user with more control over the
> implementation. Read the relevant documentation for more details.

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
manual_labeling_filepath <- "/path/to/samples_manual_labeling.csv"
samples_filepath <- "/path/to/samples.csv"
data_dir <- "/path/to/fcs_files"

single_sample_labels <- generatePopulationAssignments(
							manual_labeling_filepath, 
							samples_filepath, 
							data_dir)

# Calculate AOF for Er168Di (CD3) and Nd142Di, designating positive and negative 
# populations for each channel.

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

## CyTOF QC App

This R/Shiny application includes two features, sample background tracking and QC 
report generation, detailed in the `Features` section below. Please note: The app
is not currently hosted online; however, it can be used locally.

This app operates under the following assumptions:
- marker-containing channel names include `Di`
- DNA-containing channel names include `191Di`
- bead-contaning channel names include the following: `140Di`, `175Di`


### App Installation & Setup Instructions
1. Install R and RStudio
2. Install R package dependencies by running the following in an R session:

```r
install.packages("dplyr")
install.packages("DT")
install.packages("ggplot2")
install.packages("matrixStats")
install.packages("mclust")
install.packages("plyr")
install.packages("robustbase")
install.packages("shiny")
install.packages("shinydashboard")
install.packages("shinyFiles")
install.packages("stringr")

source("https://bioconductor.org/biocLite.R")
biocLite('flowCore')
```
3. Clone or download a copy of this repository

### App Launch Instructions
In this example, the repository was cloned to the `~/Desktop/` directory.

1. Open RStudio, change working directory to the `cytof_qc_app` location, load
the R package `shiny`, and run the application:

```r
setwd("~/Desktop/cytutils/cytof_qc_app")
library(shiny)
runApp()
```

The application will launch and you will see an interface similar to the following:

![app_launch_screenshot_img][app_launch_screenshot]

Note the two tabs on the left correspond to the two main features of the application.

### Features
#### Sample Background Tracking

This interface allows users to upload FCS files and export a separate sample 
background report per file.

The uploaded FCS files will undergo a pre-processing step that identifies beads, 
doublets, and debris. A sample background report is then generated and includes 
total event, total cell, and median channel values.

Median values for all channels present in the FCS file will be reported.

Please reference the `sample_background_report.csv` file found under 
`cytof_qc_app/docs/sample_background_report.csv` for an example report.

##### Feature Use Instructions

1. Select the directory you would like to export the sample background report to
by clicking the `Choose directory` button.

2. Upload one or multiple FCS files by clicking the `Browse...` button.

To select multiple files, hold the `Ctrl` key while clicking the files. Then click
`Open`.

A progress bar will appear on the lower-right corner.

Success or error messages will appear in the main panel in green and red text, 
respectively.

![sample_bckg_success_screenshot_img][sample_bckg_success_screenshot]

A separate csv file containing sample background information for each of the 
successfully processed files will be exported to the location you chose in step 1.

If a sample was not successfully processed, a csv file will not be exported for
the file. However, this will not affect export of sample background report files
for successfully processed files.

In addition, the the results being exported, a `Sample Background Report Outputs` 
section will appear and render the results on the screen.

#### CyTOF QC Report Generator

This interface allows users to upload FCS files and export a separate cytof 
qc report per file.

The uploaded FCS files will undergo a pre-processing step that identifies beads, 
doublets, and debris. A cytof qc report is then generated and includes 
total event, total cell, total bead, median Oxide %, and more.

Please reference the `cytof_qc_report.csv` file found under 
`cytof_qc_app/docs/cytof_qc_report.csv` for an example report.

##### Feature Use Instructions

1. Select the directory you would like to export the sample background report to
by clicking the `Choose directory` button.

2. Upload one or multiple FCS files by clicking the `Browse...` button.
![cytof_qc_file_upload_success_screenshot_img][cytof_qc_file_upload_success_screenshot]

3. Upon success, two things will happen:
	-  A `QC Report Outputs` section will render in the middle third of the app 
	screen. This will show the results that were exported.
	-  A `Gating Inspection` section will render on the lower half of the screen. 
	Select the file you would like to visualize gates for and click the `Generate 
	Gating Visualization` button.

4. If the gating visualization looks abnormal, click the `Flag Abnormal Gating`
button. Note: This action can be undone by clicking `Undo Abnormal Gating Flag`.
Doing either will update the `Abnormal Gating` column of the exported QC report.

5. You may manually adjust the gating by clicking and dragging on the gating visualization
plot. Then select the population type from the drop down (i.e. cell, bead, debris)
and click `Update Gating`.

6. When you are satisfied with how the plot looks, click the `Update QC Report` 
button to update the previously exported QC report with new values based on your 
manual gates. The table in the `QC Report Outputs` section will also be updated.

![cytof_qc_gating_screenshot_img][cytof_qc_gating_screenshot]


[app_launch_screenshot]: cytof_qc_app/docs/app_launch_screenshot.png
[sample_bckg_success_screenshot]: cytof_qc_app/docs/cytof_qc_app_sample_background_tracking.png
[cytof_qc_file_upload_success_screenshot]: cytof_qc_app_file_upload_and_report_tables.png
[cytof_qc_gating_screenshot]: cytof_qc_app/docs/cytof_qc_app_gating_inspection.png
