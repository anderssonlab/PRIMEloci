# PRIMEloci

The PRIMEloci repository focuses on genome-wide prediction of accurately identified tag clusters from CAGE data using machine learning models. The primary model was trained using light gradient boosting machine (lightGBM) on GM12878 whole-cell CAGE and nucCAGE data from the Andersson lab. This project is designed to automate the process from bigwig files to accurately identified BED files. 

Additionally, it provides flexibility for users to apply their existing .rds of RSE, GRanges, or GRangesList objects from their analyses and pass them through the model. While the project was initially developed based on the human genome hg38, it can be adapted for use with other species.

## Table of Contents

- [Overview](#overview)
- [Installation](#installation)
- [Usage](#usage)
  - [Configuration](#configuration)
  - [Running Scripts](#running-scripts)
- [Contributing](#contributing)
- [License](#license)

## Overview

PRIMEloci provides a streamlined, automated workflow for handling CAGE-seq data, focusing on accurate prediction of tag clusters using machine learning models. The process involves five key steps:

1. **Extracting CTSS Data**: Create CAGE Transcriptional Start Site (CTSS) object from bigwig files.
2. **Generating Tag Cluster Data**: Identify tag clusters (TC) from the CTSS object.
3. **Profiling TSS Clusters**: Prepare profiles from CTSS and TC data for model input.
4. **Predicting TSS Profile Probabilities**: Use PRIMEloci models to predict profile probabilities.
5. **Filtering Prediction Results**: Produce non-overlapping lists in .bed and .rds formats for further R analysis.

Scripts can be run individually or as part of a pipeline, with a main bash script for selective execution. Advanced users can also prepare data and train models to customize the workflow. The detailed steps are outlined below.

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
   install.packages("path/to/PRIMEloci_0.6.tar.gz", repos = NULL, type = "source")
   ```

   Check the DESCRIPTION file for other relevant R packages that need to be installed.

3. **Install Python Packages**:

   Ensure you have Python 3.9 or higher. Install the required Python packages:

   ```bash
   pip install matplotlib numpy pandas seaborn lightgbm scikit-learn
   ```

4. **Make the Main Script Executable**:

   ```bash
   chmod +x PRIMEloci.sh
   ```

These steps will set up the necessary environment for running PRIMEloci scripts. You can now proceed with executing the main script or individual scripts as needed.

## Usage

To use PRIMEloci, follow these steps:

1. **Navigate to the Genome-wide Prediction Directory**

   ```bash
   cd genomewide_prediction
   ```

2. **Configure Parameters**

   Modify the `config.sh` file to set all necessary variables. This file contains all the parameter settings required for the scripts to run. Ensure you set the paths for your CAGE bigWig directory and the design matrix describing the CAGE bigWig files.

3. **Run the Scripts**

   To execute the entire workflow, run:

   ```bash
   ./PRIMEloci.sh --all
   ```

   This will process the CAGE bigWig data from the initial extraction to the final output of non-overlapping lists in .bed and .rds formats.

### Examples usage in some cases

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
./PRIMEloci.sh -1 -3 -4
```

## Contributing

Contributions are welcome! Please submit a pull request or open an issue to discuss any changes or improvements.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
