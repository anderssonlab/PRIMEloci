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
conda env remove -n prime-conda-env -y >/dev/null 2>&1 || true
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
conda install -y -c conda-forge r-r.utils r-future r-future.apply r-future.callr r-foreach r-argparse r-doparallel r-reticulate r-arrow r-igraph r-catools r-zoo r-biocmanager r-remotes r-devtools
```

**Expect:** R packages installed without compilation.

---

8. System libraries & toolchain for Bioconductor builds
Some Bioconductor packages (e.g. CAGEfightR, rtracklayer) need to compile C, C++, or Fortran code, and link against system libraries (XML, compression, Unicode, etc.).
So before installing with BiocManager, make sure your environment has:
    Compilers: C, C++, Fortran
    Build tools: make, pkg-config
    Core libraries:
    libxml2, libiconv (text/XML handling)
    libcurl, openssl (network + HTTPS support)
    zlib, xz, bzip2 (compression)
    pcre2, icu (regex + Unicode support)

Example: Linux (x86_64)
```bash
conda install -y -c conda-forge gcc_linux-64 gxx_linux-64 gfortran_linux-64 make pkg-config libxml2 libiconv libcurl openssl zlib xz bzip2 pcre2 icu
# Expect inside R:
# x86_64-conda-linux-gnu-cc
# x86_64-conda-linux-gnu-c++ -std=gnu++17
# x86_64-conda-linux-gnu-gfortran
```

Example: macOS Apple Silicon (arm64)
```bash
conda install -y -c conda-forge clang_osx-arm64 clangxx_osx-arm64 gfortran_osx-arm64 llvm-openmp make pkg-config libxml2 libiconv libcurl openssl zlib xz bzip2 pcre2 icu
# If gfortran_osx-arm64 isn’t available, skip it — most packages build fine without Fortran.
# Expect inside R:
# .../envs/prime-conda-env/bin/clang
# .../envs/prime-conda-env/bin/clang++ -std=gnu++17
# .../envs/prime-conda-env/bin/gfortran
```

Example: macOS Intel (x86_64)
```bash
conda install -y -c conda-forge clang_osx-64 clangxx_osx-64 gfortran_osx-64 llvm-openmp make pkg-config libxml2 libiconv libcurl openssl zlib xz bzip2 pcre2 icu
# Expect inside R:
# .../envs/prime-conda-env/bin/clang
# .../envs/prime-conda-env/bin/clang++ -std=gnu++17
# .../envs/prime-conda-env/bin/gfortran
```

Example: Windows
On Windows, conda does not provide compilers. You need:
    Rtools (for base R + Bioconductor builds):
    Download from: https://cran.r-project.org/bin/windows/Rtools/
    Add to PATH (installer does this by default).
    System libraries are usually bundled with Rtools; if needed, install additional ones via conda-forge’s m2w64-* packages.
```bash
# Expect inside R:
# x86_64-w64-mingw32-gcc
# x86_64-w64-mingw32-g++ -std=gnu++17
# x86_64-w64-mingw32-gfortran
```

Verify compilers in R (all platforms)
```bash
R -q --vanilla -e 'Sys.getenv(c("CC","CXX","FC")); system("R CMD config CC"); system("R CMD config CXX"); system("R CMD config FC")'
```
**Expect:** Compiler paths pointing inside your conda env (Linux/macOS) or Rtools (Windows).

---

## 9. Install Bioconductor packages via conda or inside R (recommended on macOS arm64)
```bash
conda install -y -c conda-forge -c bioconda bioconductor-rtracklayer bioconductor-genomicranges bioconductor-iranges bioconductor-genomeinfodb bioconductor-summarizedexperiment bioconductor-biocparallel bioconductor-bsgenome bioconductor-cagefightr
```
or
```bash
R -q --vanilla <<'RSCRIPT'
install.packages("BiocManager", repos="https://cloud.r-project.org")
BiocManager::install(c(
  "CAGEfightR",
  "rtracklayer",
  "GenomicRanges",
  "IRanges",
  "GenomeInfoDb",
  "SummarizedExperiment",
  "BiocParallel",
  "BSgenome"
), update=TRUE, ask=FALSE)
RSCRIPT
```

**Expect:** Bioconductor resolves versions compatible with R 4.4 and installs successfully.

---

## 10. Verify Bioconductor installs


Quick check:
```bash
R -q -e 'library(CAGEfightR); packageVersion("CAGEfightR")'
```

**Expect:** Prints a version (e.g., ‘1.26.0’) with no segfaults.

---

## 9. Install compilers (for building PRIME and PRIMEloci)
On Linux, use:
```bash
conda install -y -c conda-forge gxx_linux-64 gcc_linux-64 gfortran_linux-64 make pkg-config
```

Verify inside R:
```r
R -q -e 'Sys.getenv(c("CC","CXX","FC")); system("R CMD config CC"); system("R CMD config CXX"); system("R CMD config FC")'
```

**Expect:** Compiler paths inside conda, like:
```
x86_64-conda-linux-gnu-cc
x86_64-conda-linux-gnu-c++ -std=gnu++17
x86_64-conda-linux-gnu-gfortran
```

---

## 10. Install bcp from GitHub
```bash
R -q -e 'Sys.unsetenv("R_LIBS_USER"); target_lib <- .libPaths()[1];
if (!requireNamespace("remotes", quietly=TRUE)) install.packages("remotes", lib=target_lib, repos="https://cloud.r-project.org");
remotes::install_github("swang87/bcp", lib=target_lib, upgrade="never", dependencies=TRUE, force=TRUE);
library(bcp, lib.loc=target_lib); cat("bcp loaded from: ", system.file(package="bcp"), "\n"); print(packageVersion("bcp"))'
```

**Expect:** `bcp` loads and version prints.

## 10. Install PRIME
```bash
R -q -e 'Sys.unsetenv("R_LIBS_USER"); target_lib <- .libPaths()[1];
remotes::install_github("anderssonlab/PRIME", lib=target_lib, upgrade="never", dependencies=TRUE, force=TRUE);
library(PRIME, lib.loc=target_lib); cat("PRIME loaded from: ", system.file(package="PRIME"), "\n"); print(packageVersion("PRIME"))'
```

**Expect:** `PRIME` loads and version prints.

---

---

## 12. Install PRIMEloci (local tarball)
```bash
## Change path to PRIME directiry before running this command
R -q -e 'Sys.setenv(RETICULATE_PYTHON="~/.conda/envs/prime-conda-env/bin/python3", PYTHONNOUSERSITE="1"); install.packages("/PATH/TO/PRIME/PRIMEloci_0.2.1.tar.gz", repos=NULL, type="source", lib=.libPaths()[1]); library(PRIMEloci); packageVersion("PRIMEloci"); library(reticulate); print(py_config()); cat("\nPRIMEloci loaded OK ✅\n")'
```

**Expect:**
- PRIME installs cleanly.
- `packageVersion("PRIME")` prints `0.1.1.7`.
- `py_config()` shows Python from `prime-conda-env`.
- Final message: `PRIME loaded OK ✅`.

---

## 12. Quick PRIME test
```r
R
library(GenomicRanges)
library(PRIME)
packageVersion("PRIME")
plc_focal_example <- run_PRIMEloci_focal_example(python_path = "~/.conda/envs/prime-conda-env/bin/python3")
plc_example <- run_PRIMEloci_example(python_path = "~/.conda/envs/prime-conda-env/bin/python3")
```

**Expect:** Example runs without error.

---

✅ Done — Environment `prime-conda-env` is ready.
