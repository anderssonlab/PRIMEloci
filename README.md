# PRIMEloci

The PRIMEloci repository offers tools for genome-wide prediction of accurately identified tag clusters (TCs) from CAGE data using machine learning models. The core model, based on LightGBM, was trained on GM12878 whole-cell CAGE and nucCAGE data from the Andersson lab. PRIMEloci automates the workflow from bigWig files to accurately identified TC .bed and .rds files, providing flexibility for users to either use bash scripts for a full pipeline execution or directly interact with R functions for a more programmatic approach. While the project was initially developed based on the human genome hg38, it can be adapted for use with other species.

## Table of Contents

- [Overview](#overview)
- [Installation](#installation)
- [Usage](#usage)
- [Contributing](#contributing)
- [License](#license)

## Overview

The workflow focuses on the accurate prediction of tag clusters using machine learning, encompassing five key steps:

1. **Extracting CTSS data**: Generate CAGE transcriptional start site (CTSS) objects from bigWig files.
2. **Identifying Tag Clusters**: Derive tag clusters (TCs) from the CTSS data.
3. **Processing Profiles**: Prepare processed CAGE profiles for input into the prediction model.
4. **Predicting Profile Probabilities**: Utilize PRIMEloci models to predict profile probabilities.
5. **Filtering Results**: Output non-overlapping TC lists in .bed and .rds formats for further analysis in R.

PRIMEloci offers two primary ways to utilize its genome-wide prediction:

Bash Script: Execute a full pipeline or some parts of the pipeline using pre-configured bash scripts, ideal for users who prefer a command-line interface or need to automate large-scale data processing tasks. Each step can be run individually or as part of a pipeline with a main bash script for selective execution. Advanced users can also prepare data and train models to customize the workflow [LINK to nucCAGE repo].

R Functions: Directly use R functions provided in the PRIMEloci package, allowing for more control, customization, and integration with other R-based workflows. When using the PRIMEloci R functions, users typically handle steps 1 and 2 (CTSS and TC identification using CAGEfighR functions) within R, while the run_PRIMEloci() function covers steps 3 through 5, processing existing CTSS and TC data through the prediction model and generating the final output.

## Installation

To install PRIMEloci, follow these steps:

1. **Clone the Repository**:

   ```bash
   git clone https://github.com/anderssonlab/PRIMEloci.git
   cd PRIMEloci
   ```

2. **Install R Package**:

   Ensure you have R version 4.2 or higher. Open R or RStudio and run the following command to install the PRIMEloci package locally from the provided .tar.gz file:

   ```r
   install.packages("path/to/PRIMEloci_0.7.tar.gz", repos = NULL, type = "source")
   ```

   Check the DESCRIPTION file for other relevant R packages that need to be installed.

3. **Install Python Packages**:

   Ensure you have Python 3.9 or higher. Install the required Python packages:

   ```bash
   pip3 install matplotlib numpy pandas seaborn lightgbm scikit-learn
   ```

4. **Make the Main Script Executable**:

   ```bash
   cd genomewide_prediction
   chmod +x PRIMEloci.sh
   ```

These steps will set up the necessary environment for running PRIMEloci scripts. You can now proceed with executing the main script or individual scripts as needed.

## Usage of bash script

To use PRIMEloci, follow these steps:

1. **Navigate to the Genome-wide Prediction Directory**

   ```bash
   cd genomewide_prediction
   ```

2. **Configure Parameters**

   Modify the `bash_config_PRIMEloci.sh` file to set all necessary variables. This file contains all the parameter settings required for the scripts to run. The default parameters was set to use the K562 CAGE files located in `path/to/PRIMEloci/example/resources`.

3. **Run the Scripts**

   To execute the entire workflow, run:

   ```bash
   ./PRIMEloci.sh --all
   ```

   This will process the CAGE bigWig data from the initial extraction to the final output of non-overlapping lists in .bed and .rds formats.

4. **Examples usage in some cases**

   Run all steps:

   ```bash
   ./PRIMEloci.sh --all
   ```

   Run only step as CTSS- and TC-.rds files existed:

   ```bash
   ./PRIMEloci.sh -p
   ```

   Run specific steps:

   ```bash
   ./PRIMEloci.sh -4 -5
   ```

### Usage of R function in PRIMEloci package

Ensure that all settings for run_PRIMEloci() are correctly configured in `path/to/PRIMEloci/genomewide_prediction/config_R_PRIMEloci.yaml`. 

```R
# load libraries
library(PRIMEloci)

# Load necessary R objects
ctss_rse <- loadRDS("path_to_ctss_rse.rds")
tc_grl <- loadRDS("path_to_tc_grl.rds")

# Run the PRIMEloci workflow (Steps 3-5)
primeloci_tc_grl <- run_PRIMEloci(ctss_rse = ctss_rse,
                                  tc_grl = tc_grl,
                                  config_file = "config_R_PRIMEloci.yaml")
```

The example of complete workflow, from processing bigWig files to the final steps in R, can be found in PRIMEloci_genomewide_prediction.rmd.

## Contributing

Contributions are welcome! Please submit a pull request or open an issue to discuss any changes or improvements.

## License

