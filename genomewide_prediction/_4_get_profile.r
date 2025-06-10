#!/usr/bin/env Rscript

suppressWarnings(suppressMessages({
  library(argparse)
  library(CAGEfightR)
  library(parallel)
  library(GenomicRanges)
  library(PRIME)
  library(future.apply)
  library(SummarizedExperiment)
}))

# Define argument parser
parser <- ArgumentParser()

# Input
parser$add_argument("--ctss_rse", default = "./ctss_rse.rds",
                    help = "FULL PATH to input file name for ctss_rse rds object") # nolint: line_length_linter.
parser$add_argument("--region", default = "./tc_grl.rds",
                    help = "FULL PATH to input file name for input regions rds object") # nolint: line_length_linter.

parser$add_argument("--python_path", default = "~/.virtualenvs/prime-env",
                    help = "Path to Python executable. If not provided, the system's default Python will be used.") # nolint: line_length_linter.

# Output
parser$add_argument("-o", "--output_dir", default = "./",
                    help = "Profile output directory")
parser$add_argument("--profile_dir_name", default = "PRIMEloci_profiles",
                    help = "Name of main profile directory")
parser$add_argument("-s", "--save_count_profiles", action = "store_true",
                    default = FALSE,
                    help = "Flag to save count profiles. Default is FALSE.")
parser$add_argument("-f", "--file_format", type = "character",
                    default = "npz", choices = c("parquet", "csv", "npz"),
                    help = "File format for output files. Default is 'npz'.")
parser$add_argument("--addtn_to_filename", default = "")

# Parameters
parser$add_argument("-p", "--num_cores", type = "integer", default = NULL,
                    help = "Number of cores to use for parallel processing")
parser$add_argument("-e", "--ext_dis", default = 200,
                    help = "Extension distance")

# Parse arguments
args <- parser$parse_args()

# Setup
infile_ctss_rse <- args$ctss_rse
infile_tc_grl <- args$region

output_dir <- args$output_dir
PRIME:::create_output_dir(output_dir)
PRIME::plc_setup_tmp_dir(output_dir)

profile_dir_name <- args$profile_dir_name
save_count_profiles <- args$save_count_profiles
file_format <- args$file_format

num_cores <- args$num_cores
ext_dis <- as.integer(args$ext_dis)

# Assertions
if (!is.null(num_cores)) {
  assertthat::assert_that(
    is.numeric(num_cores),
    num_cores %% 1 == 0,
    num_cores > 0,
    msg = "`num_cores` must be a positive integer or NULL."
  )
}

assertthat::assert_that(
  is.numeric(ext_dis),
  ext_dis %% 1 == 0,
  ext_dis > 0,
  msg = "`ext_dis` must be a positive integer."
)

if (is.null(num_cores)) {
  num_cores <- max(1, min(25, parallel::detectCores() %/% 2))
}
if (num_cores == 1) {
  processing_method <- "callr"
  plc_message("⚠️ num_workers was set to 1. Using callr backend: tasks will run sequentially (despite using multiple R sessions).") # nolint: line_length_linter.
} else {
  processing_method <- PRIME:::plc_detect_parallel_plan()
}

# Python config
PRIME::plc_message("🚀 Setting up Python environment")

if (is.null(args$python_path)) {
  py <- reticulate::import("sys")
  python_path <- py$executable
} else {
  python_path <- args$python_path
}
py_conf <- PRIME:::plc_configure_python(python_path = python_path)
check_npz <- PRIME::plc_test_scipy_save_npz()
if (!check_npz) {
  plc_message("⚠️ Falling back to .parquet format")
  file_type <- "parquet"
} else {
  plc_message("✅ Using .npz format")
  file_type <- "npz"
}

# Load input
writeLines("\nReading input data..")
ctss_rse <- readRDS(infile_ctss_rse)
tc_grl <- readRDS(infile_tc_grl)

# Run profiling
PRIME::plc_message("🚀 Running PRIMEloci -4: compute count & normalized profiles for each sample") # nolint: line_length_linter.
PRIME::plc_profile(
  ctss_rse,
  tc_grl,
  output_dir,
  profile_dir_name,
  file_type = file_type,
  python_path = py_conf$python,
  addtn_to_filename = args$addtn_to_filename,
  save_count_profiles = save_count_profiles,
  num_cores = num_cores,
  processing_method = processing_method,
  ext_dis
)
PRIME::plc_message("✅ DONE :: Profiles computed and saved.")
