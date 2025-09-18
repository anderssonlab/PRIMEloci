# 🚀 How to Install PRIMEloci on RStudio Server (R 4.2.2)

PRIMEloci has been developed and tested under **R 4.4** (with Bioconductor 3.19).
While it may be possible to install and run PRIMEloci under **R 4.2** (with Bioconductor 3.16), you may encounter unexpected errors due to package version differences.

## Step 1: Load Required Modules

Run the following commands in the **terminal** to load the necessary modules and open the correct version of R:

```bash
module load gcc/11.2.0
module load R/4.2.2

git clone https://github.com/anderssonlab/PRIMEloci.git

R
```

---

## Step 2: Install PRIME and PRIMEloci from GitHub

In **R**, run the following commands to install the **PRIME** package from GitHub:

```r
library(devtools)
devtools::install_github("anderssonlab/PRIME")
# Full installation protocol for PRIME can be found here:
# https://github.com/anderssonlab/PRIME/blob/main/PRIME_installation.md
```

For **PRIMEloci**, install from a local `.tar.gz` (the model will be set up in the PRIME installation directory):

```r
install.packages("/PATH/TO/PRIMEloci/PRIMEloci_1.0.tar.gz")
# Full installation protocol for PRIMEloci can be found here:
# https://github.com/anderssonlab/PRIMEloci/blob/main/INSTALL_PRIMEloci.md
```

Alternatively, if installing via devtools (then copy the model as described in step 4):

```r
devtools::install_github("anderssonlab/PRIMEloci")
```

---

## Step 3: Open RStudio Server and Setup a PRIME-Compatible Virtual Environment

1. Open RStudio Server in your browser:  
   `http://scarbrna01fl:8787/`

2. Access the PRIMEloci environment file in R:

```r
txt_file <- file.path(find.package("PRIMEloci"), "envfile", "environment.txt")
txt_content <- readLines(txt_file)

library(reticulate)

# Install Python locally
reticulate::install_python(version = "3.9.22")

# Create a new virtual environment for running PRIMEloci
virtualenv_create("PRIMEloci", python = "3.9.22")

# Activate the environment and install dependencies
use_virtualenv("PRIMEloci")
virtualenv_install(envname = "PRIMEloci", packages = txt_content)

# Optional (if not already covered by environment.txt)
virtualenv_install(envname = "PRIMEloci", packages = "matplotlib")
```

---

## Step 4: Load Required Packages and Setup the Model

Load required system libraries:

```r
dyn.load("/opt/software/netcdf-c/4.8.1/lib/libnetcdf.so.19")
dyn.load("/opt/software/libxml2/2.9.9/lib/libxml2.so")
```

Load R packages:

```r
library(PRIMEloci)
library(reticulate)
library(GenomicRanges)

# Activate virtual environment
use_virtualenv("PRIMEloci", required = TRUE)
```

Copy model to PRIMEloci model directory (skip this if installed via `.tar.gz`):

```r
current.model <- "/maps/projects/ralab/data/projects/nucleiCAGEproject/7.Model_development/PRIMEloci_GM12878_model_1.0.sav"
target.model <- file.path(find.package("PRIME"), "model", "PRIMEloci_GM12878_model_1.0.sav")

if (!file.exists(target.model)) {
  file.copy(from = current.model, to = target.model)
} else {
  message("Model file already exists at the target location. No need to copy.")
}
```

---

## Step 5: Run the Test Example

```r
ctss.file <- system.file("extdata", "ctss_rse_chr16to17.rds", package = "PRIMEloci")
ctss <- readRDS(ctss.file)
PL.test <- PRIMEloci(ctss, python_path = paste0(reticulate::virtualenv_root(), "/PRIMEloci"))
```

---

✅ If everything runs without errors, PRIMEloci is successfully installed and ready to use!
