#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(argparse)
  library(CAGEfightR)
  library(parallel)
  library(GenomicRanges)
  library(PRIME)
  library(future.apply)
  library(SummarizedExperiment)
})

# Define argument parser
parser <- ArgumentParser()

# Input
parser$add_argument("-c", "--infile_ctss_rse", default = "./ctss_rse.RDS",
                    help = "FULL PATH to input file name for ctss_rse rds object") # nolint: line_length_linter.
parser$add_argument("-t", "--infile_tc_grl", default = "./tc_grl.RDS",
                    help = "FULL PATH to input file name for tc_grl rds object")
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
parser$add_argument("-l", "--log", default = NULL,
                    help = "Log file name e.g. PRIMEloci-4.log")

# Parameters
parser$add_argument("-p", "--num_cores", type = "integer", default = NULL,
                    help = "Number of cores to use for parallel processing")
parser$add_argument("-e", "--ext_dis", default = 200,
                    help = "Extension distance")

# Parse arguments
args <- parser$parse_args()

# Setup
infile_ctss_rse <- args$infile_ctss_rse
infile_tc_grl <- args$infile_tc_grl

output_dir <- create_output_dir(args$output_dir)

#primeloci_tmp <- setup_tmp_dir(output_dir)

log <- if (is.null(args$log) || args$log == "NULL") NULL else args$log
log_target <- setup_log_target(log, output_dir)

file_format <- args$file_format
save_count_profiles <- args$save_count_profiles
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

# Python config
py_conf <- PRIME::configure_plc_python(python_path = args$python_path,
                                       log_target = log_target)

# Load input
writeLines("\nReading input data..")
ctss_rse <- readRDS(infile_ctss_rse)
tc_grl <- readRDS(infile_tc_grl)

# Logging start
plc_log("\n\n\n 🚀 Running PRIMEloci -4: compute count & normalized profiles for each sample", log_target) # nolint: line_length_linter.

# Run profiling
PRIME::plc_profile(
  ctss_rse = ctss_rse,
  tc_for_profile = tc_grl,
  outdir = output_dir,
  profile_dir_name = args$profile_dir_name,
  file_type = file_format,
  python_path = py_conf$python,
  addtn_to_filename = "",
  save_count_profiles = save_count_profiles,
  num_cores = num_cores,
  log_file = log_target,
  ext_dis = ext_dis
)

plc_log("✅ DONE :: Profiles computed and saved.", log_target)