Certainly! Here is the formatted content ready to be copied to your GitHub `README.md` file:

---

# PRIMEloci

PRIMEloci is a collection of scripts designed for processing and analyzing CAGE-seq data, predicting transcription start site (TSS) profiles, and filtering prediction results. This repository provides a flexible and modular approach to run various steps in the analysis pipeline, allowing users to execute specific steps or the entire pipeline based on their needs.

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

PRIMEloci provides a set of scripts for:

1. Extracting CAGE-seq data.
2. Generating TSS cluster data.
3. Profiling TSS clusters.
4. Predicting TSS profile probabilities using a machine learning model.
5. Filtering prediction results.

The scripts are designed to be run individually or as part of a pipeline, with flexibility provided via a main bash script that allows for selective execution of each step.

## Installation

Clone the repository:

```bash
git clone https://github.com/yourusername/PRIMEloci.git
cd PRIMEloci
```

Make the main script executable:

```bash
chmod +x run_scripts.sh
```

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
