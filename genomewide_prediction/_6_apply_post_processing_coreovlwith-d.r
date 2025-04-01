#!/usr/bin/env Rscript

writeLines("\n### Running coreovlwith-d postprocess ###")
writeLines("\n# Importing R libraries..")

suppressPackageStartupMessages({
  library(argparse)
  library(GenomicRanges)
  library(parallel)
  library(tools)
  library(stringr)
  library(PRIMEloci)
  library(dplyr)
  library(data.table)
  library(IRanges)
})

# Argument parser setup
parser <- ArgumentParser()

parser$add_argument("-i", "--input_file", required = TRUE,
                    help = "Path to the input BED file.")
parser$add_argument("-t", "--score_threshold", type = "double", default = 0.75,
                    help = "Score threshold for filtering GRanges.")
parser$add_argument("-o", "--output_dir", default = NULL,
                    help = "Directory to save the output files.")
parser$add_argument("-m", "--use_max_cores", action = "store_true",
                    help = "Flag to use the maximum number of available cores.")
parser$add_argument("-d", "--score_diff", type = "double", default = 0.10,
                    help = "Score difference threshold for merging.")
parser$add_argument("--partial_name", default = "pred_all",
                    help = "Partial name to search for in BED files.")

args <- parser$parse_args()

# Set parameters
score_threshold <- args$score_threshold
score_diff <- args$score_diff
bed_file <- args$input_file
use_max_cores <- args$use_max_cores
output_dir <- args$output_dir
core_width <- 151

#Check numeric parameters
assertthat::assert_that(
  is.numeric(score_threshold),
  score_threshold > 0,
  score_threshold < 1,
  msg = "`score_threshold` must be a numeric value between 0 and 1."
)
assertthat::assert_that(
  is.numeric(score_diff),
  score_diff >= 0,
  score_diff < score_threshold,
  msg = "`score_diff` must be a non-negative numeric value and smaller than `score_threshold`." # nolint: line_length_linter.
)
if (!is.null(num_cores)) {
  assertthat::assert_that(
    is.numeric(num_cores),
    num_cores %% 1 == 0,
    num_cores > 0,
    msg = "`num_cores` must be a positive integer or NULL."
  )
}
output_dir <- create_output_dir(args$output_dir)
primeloci_tmp <- setup_tmp_dir(output_dir)

log <- if (args$log == "NULL") NULL else args$log
log_target <- setup_log_target(log, output_dir)

postprocess_partial_name <- args$partial_name

# Begin postprocessing
plc_log("\n\n\n 🚀 Running PRIMEloci -6: Postprocessing prediction BEDs", log_target)

bed_files <- find_bed_files_by_partial_name(
  primeloci_tmp,
  partial_name = postprocess_partial_name,
  log_file = log_target
)

if (length(bed_files) == 0) {
  msg <- paste("❌ No BED files found for postprocessing in", primeloci_tmp)
  plc_log(msg, log_target, level = "❌ ERROR")
  stop(msg)
}

plc_log(sprintf("📂 Found %d BED file(s) for processing.", length(bed_files)), log_target)

lapply(seq_along(bed_files), function(i) {
  bed_file <- bed_files[i]
  basename_raw <- tools::file_path_sans_ext(basename(bed_file))
  pattern_match <- sub(
    paste0("^.*", postprocess_partial_name, "_(.*?)_combined.*$"),
    "\\1",
    basename_raw
  )

  sample_name <- if (identical(pattern_match, basename_raw)) {
    basename_raw
  } else {
    pattern_match
  }

  result_gr <- coreovl_with_d(
    bed_file = bed_file,
    score_threshold = score_threshold,
    score_diff = score_diff,
    core_width = core_width,
    return_gr = FALSE,
    output_dir = output_dir,
    num_cores = num_cores,
    log_file = log_target
  )

  if (!is.null(result_gr)) {
    list(name = sample_name, gr = result_gr)
  } else {
    plc_log(paste("⚠️ Skipped due to failure:", bed_file), log_target)
    NULL
  }
})

plc_log("✅ DONE :: Postprocessing completed.", log_target)
