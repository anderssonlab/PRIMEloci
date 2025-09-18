# PRIMEloci Environment Setup Guide

This guide sets up PRIME in a clean, reproducible `conda` environment.

```bash
git clone https://github.com/anderssonlab/PRIMEloci.git
```

---



## 1. Clean start

```bash
module purge
module load miniconda/24.5.0
conda config --set channel_priority strict
```

**Expect:** No errors. `conda config` confirms channel priority is set.

---



## 2. Remove old env (if exists) and create new one

```bash
# remove old env (if exists)
conda env remove -n prime-conda-env -y >/dev/null 2>&1 || true

# create the new env
conda create -y -n prime-conda-env -c conda-forge python=3.11 r-base=4.4
```

**Expect:** Conda solves environment and installs Python 3.11.x + R 4.4.x.

---



## 3. Activate & sanitize

```bash
conda activate prime-conda-env
unset LD_PRELOAD
unset LD_LIBRARY_PATH
```

**Expect:** Prompt changes to `(prime-conda-env)`.

---



## 4. Ignore `~/.local` Python site-packages inside this env

```bash
conda env config vars set PYTHONNOUSERSITE=1
conda deactivate && conda activate prime-conda-env
```

**Expect:** Env reactivates cleanly. No `.local` packages will be seen.

---



## 5. Make R prefer this env’s library automatically

```bash
RPROFILE_SITE="$CONDA_PREFIX/lib/R/etc/Rprofile.site"
mkdir -p "$(dirname "$RPROFILE_SITE")"
cat > "$RPROFILE_SITE" <<'RPROF'
## Prefer the conda env's library first
conda_lib <- "~/.conda/envs/prime-conda-env/lib/R/library"
if (dir.exists(path.expand(conda_lib))) {
  .libPaths(unique(c(path.expand(conda_lib), .libPaths())))
}
## Avoid user lib overriding this env
Sys.unsetenv("R_LIBS_USER")
RPROF
```

Verify:
```bash
R -q -e "print(.libPaths())"
```

**Expect:** First entry is `~/.conda/envs/prime-conda-env/lib/R/library`.

---



## 6. Install Python dependencies

```bash
## Change path to PRIMEloci directiry before running this command
conda install -y -c conda-forge --update-specs --file /PATH/TO/PRIME/inst/envfile/environment.txt
```

**Expect:** Conda installs specified versions based on the environment file.

Quick check:
```bash
python - <<'PY'
import sys, numpy, scipy, pandas, sklearn, joblib, pyarrow, fastparquet, lightgbm
print("python", sys.version.split()[0])
print("numpy", numpy.__version__)
print("scipy", scipy.__version__)
print("pandas", pandas.__version__)
print("scikit-learn", sklearn.__version__)
print("joblib", joblib.__version__)
print("pyarrow", pyarrow.__version__)
print("fastparquet", fastparquet.__version__)
print("lightgbm", lightgbm.__version__)
PY
```

---



## 7. Install prebuilt CRAN R packages

```bash
conda install -y -c conda-forge r-r.utils r-future r-future.apply r-future.callr r-foreach r-argparse r-doparallel r-reticulate r-arrow r-igraph r-catools r-zoo r-biocmanager r-remotes r-devtools r-sparsematrixstats
```

**Expect:** R packages installed without compilation.

---



## 8. Install compilers

Some Bioconductor packages (e.g. **CAGEfightR**, **rtracklayer**) need to compile C, C++, or Fortran code, and link against system libraries (XML, compression, Unicode, etc.).

Follow the instructions for your operating system and R setup.

✅ **Linux (conda-based R)**

```bash
conda install -y -c conda-forge gxx_linux-64 gcc_linux-64 gfortran_linux-64 make pkg-config
```

✅ **macOS (conda-based R)**

If you run **r-base from conda** (not CRAN/Homebrew R), use the macOS-specific compilers from conda-forge:

​	**Apple Silicon (arm64):**

```
conda install -y -c conda-forge \
  clang_osx-arm64 clangxx_osx-arm64 gfortran_osx-arm64 llvm-openmp \
  make pkg-config libxml2 libiconv libcurl openssl zlib xz bzip2 pcre2 icu
```

**	Intel (x86_64):**

```
conda install -y -c conda-forge \
  clang_osx-64 clangxx_osx-64 gfortran_osx-64 llvm-openmp \
  make pkg-config libxml2 libiconv libcurl openssl zlib xz bzip2 pcre2 icu
```

🔎 **Verify compilers inside R**

```bash
R -q -e 'Sys.getenv(c("CC","CXX","FC")); 
         system("R CMD config CC"); 
         system("R CMD config CXX"); 
         system("R CMD config FC")'
```

**Expect:**

- Linux + conda → `x86_64-conda-linux-gnu-*` compilers
- macOS + conda R → `clang_osx-*` compilers

---



## 9. Install Bioconductor packages (Conda binaries)

```bash
conda install -y -c conda-forge -c bioconda bioconductor-rtracklayer bioconductor-genomicranges bioconductor-iranges bioconductor-genomeinfodb bioconductor-summarizedexperiment bioconductor-biocparallel bioconductor-bsgenome bioconductor-cagefightr
```

Quick check:

```bash
R -q -e 'library(CAGEfightR); packageVersion("CAGEfightR")'
```

**Expect:** Prints a version (e.g., ‘1.26.0’) with no segfaults.

---



## 10. Install bcp from GitHub

```bash
R -q -e 'Sys.unsetenv("R_LIBS_USER"); target_lib <- .libPaths()[1];
if (!requireNamespace("devtools", quietly=TRUE)) install.packages("devtools", lib=target_lib, repos="https://cloud.r-project.org");
devtools::install_github("swang87/bcp", lib=target_lib, upgrade="never", dependencies=TRUE, force=TRUE);
library(bcp, lib.loc=target_lib); cat("bcp loaded from: ", system.file(package="bcp"), "\n"); print(packageVersion("bcp"))'
```

**Expect:** `bcp` loads and version prints.

---



## 10. Install PRIME from GitHub

```bash
R -q -e 'Sys.unsetenv("R_LIBS_USER"); target_lib <- .libPaths()[1];
if (!requireNamespace("devtools", quietly=TRUE)) install.packages("devtools", lib=target_lib, repos="https://cloud.r-project.org");
devtools::install_github("anderssonlab/PRIME", lib=target_lib, upgrade="never", dependencies=TRUE, force=TRUE);
library(PRIME, lib.loc=target_lib); cat("PRIME loaded from: ", system.file(package="PRIME"), "\n"); print(packageVersion("PRIME"))'
```

**Expect:** `PRIME` loads and version prints.

---



## 11. Install PRIMEloci from local tarball

```bash
R -q -e 'Sys.setenv(RETICULATE_PYTHON="~/.conda/envs/prime-conda-env/bin/python3", PYTHONNOUSERSITE="1"); install.packages("/PATH/TO/PRIMEloci/PRIMEloci_1.0.tar.gz", repos=NULL, type="source", lib=.libPaths()[1]); library(PRIMEloci); packageVersion("PRIMEloci"); library(reticulate); print(py_config()); cat("\nPRIMEloci loaded OK ✅\n")'
```

**Expect:**

- PRIME installs cleanly.
- `packageVersion("PRIMEloci")` prints `1.0`.
- `py_config()` shows Python from `prime-conda-env`.
- Final message: `PRIMEloci loaded OK ✅`.

---



## 12. Quick PRIMEloci test

```r
R

library(GenomicRanges)
library(PRIMEloci)
packageVersion("PRIMEloci")

plc_focal_example <- run_PRIMEloci_focal_example(python_path = "~/.conda/envs/prime-conda-env/bin/python3")
plc_example <- run_PRIMEloci_example(python_path = "~/.conda/envs/prime-conda-env/bin/python3")
```

**Expect:** Example runs without error.

---

✅ Done — Environment `prime-conda-env` is ready.
