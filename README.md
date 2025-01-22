
# PRIMEloci

The PRIMEloci repository offers tools for genome-wide prediction of regulatory elements from CAGE data using machine learning models. The core model, based on LightGBM, was trained on GM12878 whole-cell CAGE and nucCAGE data from the Andersson lab. PRIMEloci automates the workflow from bigWig files to accurately identified enhancers and promoters, providing flexibility for users to either use bash scripts for a full pipeline execution or directly interact with R functions for a more programmatic approach. While the project was initially developed based on the human genome hg38, it can be adapted for use with other species.



## Table of Contents

- [Overview](#overview)
- [Installation](#installation)
- [Usage](#usage)
- [Contributing](#contributing)
- [License](#license)



## Overview

The workflow focuses on the prediction of regulatory elements using machine learning, encompassing six key steps:

1. **Extracting CTSS Data**: Extract CAGE transcriptional start site (CTSS) objects from bigWig files using the CAGEfightR package.
2. **Identifying Tag Clusters (TCs)**: Identify tag clusters (TCs) from the extracted CTSS data using the CAGEfightR package.
3. **Sliding Through TCs**: Slide through the identified TCs, default setting the window size to 20, to prepare data for downstream analysis.
4. **Creating Normalized Profiles**: Generate normalized profiles for input into the prediction model.
5. **Predicting Profile Probabilities**: Use PRIMEloci models to predict the probabilities of regulatory elements.
6. **Post-Processing**: Refine and filter predictions using additional criteria for improved accuracy, outputting non-overlapping loci in `.bed` format for further analysis in R.


PRIMEloci offers two primary ways to utilize its genome-wide prediction:

### Bash Script

Execute a full pipeline or some parts of the pipeline using pre-configured bash scripts, ideal for users who prefer a command-line interface or need to automate large-scale data processing tasks. Each step can be run individually or as part of a pipeline with a main bash script for selective execution. Advanced users can also prepare data and train models to customize the workflow.

### R Functions

Directly use R functions provided in the PRIME package, allowing for more control, customization, and integration with other R-based workflows. When using the PRIMEloci functions, users typically handle steps 1 and 2 (CTSS and TC identification using CAGEfighR functions) within R, while the `PRIME::PRIMEloci()` function covers steps 3 through 6, processing existing CTSS and TC data through the prediction model and generating the final output.



## Installation

To install PRIMEloci, follow these steps:

1. **Clone the Repository**:

   ```bash
   git clone git@github.com:anderssonlab/PRIMEloci.git
   cd PRIMEloci
   ```

2. **Install R Package**:

   Ensure you have R version 4.4 or higher. Open R or RStudio and run the following command to install the PRIMEloci package locally from the provided `.tar.gz` file:

   ```bash
   # KUIT-server user:
   module load openjdk/20.0.0 gcc/13.2.0 R/4.4.0
   module load hdf5 netcdf-c/4.8.1

   # Access R
   R
   ```

   ```r
   install.packages("PRIMEloci_0.9.7.tar.gz", repos = NULL, type = "source", lib="/local/path")
   ```

3. **Install Python Packages**:

   Ensure you have Python 3.9 or higher. Install the required Python packages:

   ```bash
   # KUIT-server user:
   module load python/3.9.9 python-packages/3.9

   # Install Python Packages
   pip3 install matplotlib numpy pandas seaborn lightgbm scikit-learn
   ```

4. **Make the Main Script Executable**:

   ```bash
   cd genomewide_prediction
   chmod +x PRIMEloci.sh
   ```
5. **Copy the example resources and model for demonstration**:

   ```bash
   # KUIT-server user:
   cd ../example/resources/cage_bw
   cp /maps/projects/ralab/data/projects/nucleiCAGEproject/8.Genomewide_prediction/example_cage_K562_bw/* .

   cd ../../../model
   cp /maps/projects/ralab/data/projects/nucleiCAGEproject/7.Model_development/PRIMEloci_GM12878_model_1.0.sav . 

   ``` 

These steps will set up the necessary environment for running PRIMEloci scripts. You can now proceed with executing the main script or individual scripts as needed.



## Usage of Bash Script

To use PRIMEloci, follow these steps:

1. **Navigate to the Genome-wide Prediction Directory**

   ```bash
   cd PRIMEloci/genomewide_prediction
   ```

2. **Configure Parameters**

   The example of config file can be found at `path/to/PRIMEloci/genomewide_prediction/bash_config_PRIMEloci.sh` Provide a valid configuration file with `--config`. This file must define all required parameters for the script. Example settings support CAGE files located in `PRIMEloci/example/resources`.

3. **Run the Scripts**

   Execute the desired steps with the following options:

   - **Run All Steps**: Process data from initial CTSS extraction to final post-processing of regulatory element predictions.
     ```bash
     ./PRIMEloci.sh --config <config_file> --all

     # example <config_file>
     ./PRIMEloci.sh --config bash_config_PRIMEloci.sh --all
     ```
     If server storage is not a concern, it is recommended to use --keep_tmp, as it allows you to retain temporary files for further analysis.
     ```bash
     ./PRIMEloci.sh --config <config_file> --all --keeptmp
     ```    

   - **Pre-processed Data**: Skip earlier steps if `.rds` files for CTSS and TCs/regions are already available. 
     Make sure that the files are in the correct format. Examples can be explored by running the provided example commands with `--all`. 
     ```bash
     ./PRIMEloci.sh --config <config_file> --pred  
     ```
     If server storage is not a concern, it is recommended to use --keep_tmp, as it allows you to retain temporary files for further analysis.
     ```bash
     ./PRIMEloci.sh --config <config_file> --pred --keeptmp
     ```

   - **Specific Steps**: Run individual steps as needed. 
   Note that when running step by step, all temporary files will be kept, allowing you to inspect them. You can remove these files later from
     ```bash
     ./PRIMEloci.sh --config <config_file> -3 -4 
     ```
   Each step corresponds to a specific script or function:
   - `-1`: Extract CTSS data.
   - `-2`: Identify tag clusters.
   - `-3`: Slide through TCs with a specific window size.
   - `-4`: Create normalized profiles.
   - `-5`: Predict regulatory element probabilities.
   - `-6`: Apply post-processing to refine predictions.

### Additional Applications of PRIMEloci

[will be updated]



### Usage of R Function in PRIMEloci Package

[will be changed] Ensure that all settings for `PRIMEloci()` are correctly configured in `path/to/PRIMEloci/genomewide_prediction/config_R_PRIMEloci.yaml`. 

```R
# load libraries
library(PRIMEloci)

# Load necessary R objects
ctss_rse <- loadRDS("path/to/ctss_rse.rds")
tc_grl <- loadRDS("path/to/tc_grl.rds")

# Run the PRIMEloci workflow (Steps 3-5)
primeloci_tc_grl <- run_PRIMEloci(ctss_rse = ctss_rse,
                                  tc_grl = tc_grl,
                                  config_file = "config_R_PRIMEloci.yaml")
```

## Contributing

Contributions are welcome! Please submit a pull request or open an issue to discuss any changes or improvements.

## License
