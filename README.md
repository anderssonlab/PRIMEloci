# PRIMEloci

The PRIMEloci repository focuses on genome-wide prediction of accurately identified tag clusters from CAGE data using machine learning models. The primary model was trained using light gradient boosting machine (lightgbm) on GM12878 whole-cell CAGE and nucCAGE data from the Andersson lab. This project is designed to automate the process from bigwig files to accurately identified BED files. Additionally, it provides flexibility for users to apply their existing .rds GRanges (GRangesList) objects from their analyses and pass them through the model.

## Table of Contents

- [Overview](#overview)
- [Installation](#installation)
- [Usage](#usage)
  - [Configuration](#configuration)
  - [Running Scripts](#running-scripts)
- [Scripts](#scripts)
  - [get_ctss_from_bw.R](#1-get_ctss_from_bwr)
  - [get_tc_grl.R](#2-get_tc_grlr)
  - [get_tc_profiles.R](#3-get_tc_profilesr)
  - [predict_profile_probabilities.py](#4-predict_profile_probabilitiespy)
  - [filter_bed_to_reduced_gr.R](#5-filter_bed_to_reduced_grr)
- [Contributing](#contributing)
- [License](#license)

## Overview

PRIMEloci provides a streamlined, automated workflow for handling CAGE-seq data, focusing on accurate prediction of tag clusters using machine learning models. The process involves five key steps:

1. **Extracting CAGE-seq Data**: Create TSS objects from bigwig files.
2. **Generating TSS Cluster Data**: Identify tag clusters (TC) from the TSS objects.
3. **Profiling TSS Clusters**: Prepare profiles from TSS and TC data for model input.
4. **Predicting TSS Profile Probabilities**: Use machine learning models to predict TSS profile probabilities.
5. **Filtering Prediction Results**: Produce non-overlapping lists in .bed and .rds formats for further R analysis.

Scripts can be run individually or as part of a pipeline, with a main bash script for selective execution. Advanced users can also prepare data and train models to customize the workflow. The detailed steps are outlined below.

## Installation

To install PRIMEloci, follow these steps:

1. **Clone the Repository**:

   ```bash
   git clone https://github.com/yourusername/PRIMEloci.git
   cd PRIMEloci
   ```

2. **Install R Package**:

   Ensure you have R version 4.2 or higher. Open R or RStudio and run the following command to install the PRIMEloci package locally from the provided .tar.gz file:

   ```r
   install.packages("path/to/PRIMEloci_x.x.x.tar.gz", repos = NULL, type = "source")
   ```

   Check the DESCRIPTION file for other relevant R packages that need to be installed.

3. **Install Python Packages**:

   Ensure you have Python 3.9 or higher. Install the required Python packages:

   ```bash
   pip install matplotlib numpy pandas seaborn lightgbm scikit-learn
   ```

4. **Make the Main Script Executable**:

   ```bash
   chmod +x run_scripts.sh
   ```

These steps will set up the necessary environment for running PRIMEloci scripts. You can now proceed with executing the main script or individual scripts as needed.

## Usage

### Configuration

All the parameter settings are stored in the `config.sh` file. Modify this file to match your data paths and desired parameters.

### Running Scripts

To run the scripts, use the `run_scripts.sh` file. This script accepts options to specify which steps to run.

#### Examples

Run all steps:

```bash
./run_scripts.sh --all
```

Run specific steps:

```bash
./run_scripts.sh -1 -3 -4
```

### Script Details

#### 1. get_ctss_from_bw.R

Extracts CAGE-seq data from bigWig files.

**Usage:**

```bash
Rscript _get_ctss_from_bw.r -i <CAGE_DIR> -m <DESIGN_MATRIX> -o <OUTPUT_DIR> -c <CTSS_RSE_NAME> -k
```

#### 2. get_tc_grl.R

Generates TSS cluster data from extracted CAGE-seq data.

**Usage:**

```bash
Rscript _get_tc_from_ctss.r -c <OUTPUT_DIR>/<CTSS_RSE_NAME> -o <OUTPUT_DIR> -t <TC_GRL_NAME> -e <EXTENSION_DISTANCE>
```

#### 3. get_tc_profiles.R

Profiles the TSS clusters.

**Usage:**

```bash
Rscript _get_tc_profiles.r -c <OUTPUT_DIR>/<CTSS_RSE_NAME> -t <OUTPUT_DIR>/<TC_GRL_NAME> -o <OUTPUT_DIR> -n <PROFILE_MAIN_DIR> -r <PROFILE_SUB_DIR>
```

#### 4. predict_profile_probabilities.py

Predicts TSS profile probabilities using a pre-trained model.

**Usage:**

```bash
python3 _predict_profile_probabilities.py -w <SCRIPT_DIR> -m <MODEL_PATH> -p <OUTPUT_DIR>/<PROFILE_MAIN_DIR> -r <PROFILE_SUB_DIR> -n <PREFIX_OUT_NAME> -t <THRESHOLD>
```

#### 5. filter_bed_to_reduced_gr.R

Filters prediction results to produce a reduced set of genomic ranges.

**Usage:**

```bash
Rscript _filter_bed_to_reduced_gr.r -i <FILE>
```

## Contributing

Contributions are welcome! Please submit a pull request or open an issue to discuss any changes or improvements.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

You can copy and paste this directly into your `README.md` file on GitHub. If you need any additional modifications or have further questions, feel free to ask!
