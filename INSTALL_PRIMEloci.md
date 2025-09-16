# PRIME Installation Preparation Guide

This guide explains how to set up the required **R** and **Python** environments to use the `PRIMEloci` R package.

```bash
# In terminal
git clone https://github.com/anderssonlab/PRIMEloci.git
```

For macOS users: libomp is required for LightGBM to enable OpenMP (multithreading). Without libomp, LightGBM may fail to use multithreading properly and can produce silent errors or crashes during training without clear messages. (To follow this setup, Xcode Command Line Tools and Homebrew are required.)

```bash
# Check if libomp exists

## Apple Silicon (M1/M2/M3/...)
ls /opt/homebrew/opt/libomp/lib/libomp.dylib
## Intel macOS (x86_64)
ls /usr/local/opt/libomp/lib/libomp.dylib
```

```bash
# If libomp does not exist

brew install libomp
```
On macOS, R (via homebrew libomp) and Python environments (via conda or virtualenv libomp) can conflict when used together through reticulate, potentially causing crashes with errors. If this occurs, managing environment variables or aligning libomp paths may be necessary to avoid conflicts while maintaining multithreaded performance.

---

## R Installation for PRIMEloci

The R package `PRIMEocil` depends on a mix of CRAN and Bioconductor packages, and a custom GitHub version of `PRIME`.

Note that we recommend using R 4.4 or higher. While R versions >4.2 can also be used, they may require additional setup steps. For example, on macOS, you may need to run:
```bash
# core build tools for R packages
brew install gcc pkg-config

# libraries for graphics, fonts, and rendering
brew install freetype harfbuzz fribidi libpng cairo

# libraries for image support
brew install jpeg libpng libtiff

# version control and Git support
brew install libgit2
```
when using R 4.2 to install system libraries required for building certain R packages from source.

Optional (Recommended for macOS users on R 4.2.x):
To avoid X11-related warnings and enable x11() graphics, install XQuartz. This is not required if you only use RStudio or file-based plots, but it ensures maximum compatibility with all packages and plotting functions in R.

### ✅ Full R Setup

```bash
R
```

```r
# 1. Install required CRAN packages
install.packages(c(
  "R.utils",
  "assertthat",
  "data.table",
  "future",
  "future.apply",
  "future.callr",
  "foreach",
  "argparse",
  "doParallel",
  "reticulate",
  "arrow",
  "stringr",
  "parallel",
  "magrittr"
  "tools"
))

# 2. Install BiocManager (if not already installed)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# 3. Install required Bioconductor packages
## General principle: Installing `CAGEfightR` with:

BiocManager::install("CAGEfightR")

## will automatically pull in dependencies.
## You do NOT need to install these manually unless an error occurs.**

## If errors occur, install in layers:
## core:
BiocManager::install("S4Vectors", "IRanges", "GenomeInfoDb")
## data structures:
BiocManager::install("SummarizedExperiment", "GenomicRanges", "sparseMatrixStats")
## utilities
BiocManager::install("BiocParallel", "BSgenome", "rtracklayer")
## CAGEfightR
BiocManager::install("CAGEfightR")

# 4. Install devtools if needed
if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")

# 6. Install PRIME
## Install from .tar.gz (model will be set up at PRIME inst directory),
install.packages("/PATH/TO/PRIME/PRIME_0.1.1.6.tar.gz")

# 7. Install PRIMEloci
## Install from .tar.gz (model will be set up at PRIME inst directory),
install.packages("/PATH/TO/PRIME/PRIME_0.1.1.6.tar.gz")

## otherwise, make sure that the model in the path was exist.
## The published PRIMEloci model can be found at XXXXXXXXXXXXXXXXXXXXX.
```

---

## 🐍 Python Environment Setup (Compatible with PRIMEloci functions)

You can prepare the Python environment in **four ways**. The recommended default is the `reticulate`-managed virtualenv using **NumPy 2.x** (use NumPy 1.x only if required for legacy compatibility)

---

### 🌐 Option 1 : Virtualenv via `reticulate` in R

This is the **default and recommended** method for setting up the Python environment using only R. It ensures full compatibility with `PRIME` and `PRIMEloci()` 

Reticulate virtualenv (in R) is easiest for R-focused workflows without needing external software, **but requires running `use_virtualenv()` in each session** and is not ideal for CLI use.

#### Setup instructions in R:

```bash
cd PRIME
R
```

```r
library(reticulate)

# Define the environment name
env_name <- "prime-env"

# Create the virtual environment if it doesn't exist
virtualenv_create(envname = env_name)

# Use and configure it
use_virtualenv(env_name, required = TRUE)

# Install Python packages
required_pkgs <- readLines("/PATH/TO/PRIME/inst/envfile/environment.txt")
reticulate::py_install(packages = required_pkgs, envname = env_name, method = "virtualenv")

# Confirm active Python path
py_config()$python
```

#### Example usage:

```r
library(PRIME)
library(GenomicRanges)
run_PRIMEloci_focal_example(python_path = py_config()$python)
run_PRIMEloci_example(python_path = py_config()$python)
```

---


### 🔧 Option 4: Use existing Python installation with `requirements.txt`

If you already have a working Python installation and want to use it directly:

```bash
pip3 install --upgrade pip
pip3 install -r inst/envfile/environment.txt
which python3

# Use this path in your R function call
```

```r
library(PRIME)
library(GenomicRanges)
run_PRIMEloci_focal_example(python_path = "/path/to/your/python")
run_PRIMEloci_example(python_path = "/path/to/your/python")
```

### 🌬️ Option 2: Conda

`environment.yml` is included in the `inst/envfile` folder.

#### How to create the environment:

```bash
cd PRIME

# Recommended: NumPy 2.x
conda env create -f inst/envfile/environment.yml

# Activate the environment
conda activate prime-env
which python3
conda deactivate

# Copy this path to use as python_path in R when calling run_PRIMEloci_example() or run_PRIMEloci_facet_example()
```

#### Example in R:

```r
library(PRIME)
library(GenomicRanges)
run_PRIMEloci_focal_example(python_path = "~/.conda/envs/prime-env/bin/python3")
run_PRIMEloci_example(python_path = "~/.conda/envs/prime-env/bin/python3")
```

---

### 🧪 Option 3: Virtualenv (manual setup)

This is an advanced manual option for setting up the environment outside of R. It is useful if managing environments via shell or external tools.

#### Setup instructions

```bash
python3 -m venv ~/prime_env
source ~/prime_env/bin/activate
pip3 install --upgrade pip
pip3 install -r inst/envfile/environment.txt
which python3

# Copy this path to use as python_path in R when calling run_PRIMEloci_example() or run_PRIMEloci_facet_example()
```

#### Example in R:

```r
library(PRIME)
library(GenomicRanges)
run_PRIMEloci_focal_example(python_path = "~/prime_env/bin/python3")
run_PRIMEloci_example(python_path = "~/prime_env/bin/python3")
```

---


---

## 🧠 Tips

- Always check the current Python path with `which python3` **after** activating your environment.
- Use that full path in the `python_path` argument in R.
- All four environments (reticulate virtualenv, conda, manual virtualenv, and existing Python) are compatible with `PRIMEloci()`.

---

© 2025 PRIME setup protocol | Supports PRIMEloci-compatible Python integration
