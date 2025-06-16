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

# Argument parser setup
parser <- ArgumentParser()

parser$add_argument("-i", "--input_dir", default = "./",
                    help = "Path to the input BED file.")
parser$add_argument("--partial_name", default = "pred_all",
                    help = "Partial name to search for in BED files.")
parser$add_argument("-o", "--output_dir", default = "./",
                    help = "Directory to save the output files.")

parser$add_argument("-t", "--score_threshold", type = "double", default = 0.75,
                    help = "Score threshold for filtering GRanges.")
parser$add_argument("-d", "--score_diff", type = "double", default = 0.10,
                    help = "Score difference threshold for merging.")
parser$add_argument("--core_width", type = "integer", default = 151,
                    help = "Width of the core region to consider for overlaps.")

parser$add_argument("-p", "--num_cores", type = "integer", default = NULL,
                    help = "Number of cores to use for parallel processing")


args <- parser$parse_args()

# Set parameters

input_dir <- args$input_dir
if (!dir.exists(input_dir)) {
  PRIME::plc_error("❌ Input directory does not exist.")
}

postprocess_partial_name <- args$partial_name
output_dir <- args$output_dir
PRIME::plc_create_output_dir(output_dir)

score_threshold <- args$score_threshold
score_diff <- args$score_diff
core_width <- args$core_width

#Check numeric parameters
assertthat::assert_that(
  is.numeric(score_threshold),
  score_threshold > 0,
  score_threshold < 1,
  msg = "❌ `score_threshold` must be a numeric value between 0 and 1."
)

assertthat::assert_that(
  is.numeric(score_diff),
  score_diff >= 0,
  score_diff < score_threshold,
  msg = "❌ `score_diff` must be a non-negative numeric value and smaller than `score_threshold`." # nolint: line_length_linter.
)

num_cores <- args$num_cores

if (!is.null(num_cores)) {
  assertthat::assert_that(
    is.numeric(num_cores),
    num_cores %% 1 == 0,
    num_cores > 0,
    msg = "❌ `num_cores` must be a positive integer or NULL."
  )
}

if (is.null(num_cores)) {
  num_cores <- max(1, min(25, parallel::detectCores() %/% 2))
}
if (num_cores == 1) {
  processing_method <- "callr"
  PRIME::plc_message("⚠️ num_workers was set to 1. Using callr backend: tasks will run sequentially (despite using multiple R sessions).") # nolint: line_length_linter.
} else {
  processing_method <- PRIME::plc_detect_parallel_plan()
}


plc_message("🚀 Running PRIMEloci: Postprocessing prediction BEDs")
bed_files <- PRIME::plc_find_bed_files_by_partial_name(input_dir,
                                                       partial_name = postprocess_partial_name) # nolint: line_length_linter.
if (length(bed_files) == 0) {
  plc_error(paste("❌ No BED files found for postprocessing in",
                  input_dir))
}
plc_message(sprintf("📂 Found %d BED file(s) for processing.",
                    length(bed_files)))

result_named_list <- lapply(seq_along(bed_files), function(i) {
  bed_file <- bed_files[i]
  basename_raw <- tools::file_path_sans_ext(basename(bed_file))
  pattern_match <- sub(paste0("^.*",
                              postprocess_partial_name,
                              "_(.*?)_combined.*$"),
                       "\\1",
                       basename_raw)
  sample_name <- if (identical(pattern_match, basename_raw)) {
    basename_raw
  } else {
    pattern_match
  }
  result_gr <- PRIME::plc_coreovl_with_d(bed_file = bed_file,
                                         score_threshold = score_threshold,
                                         score_diff = score_diff,
                                         core_width = core_width,
                                         return_gr = TRUE,
                                         output_dir = output_dir,
                                         save_rds = TRUE,
                                         num_cores = num_cores,
                                         processing_method = processing_method)
  if (!is.null(result_gr)) {
    list(name = sample_name, gr = result_gr)
  } else {
    plc_message(paste("⚠️ Skipped due to failure:", bed_file))
    NULL
  }
})

# Filter out failed/null entries
result_named_list <- Filter(Negate(is.null), result_named_list)
if (length(result_named_list) == 0) {
  plc_error("❌ All postprocessing attempts failed or returned NULL.")
}
plc_message(sprintf("✅ DONE :: Postprocessed %d file(s) successfully.",
                    length(result_named_list)))
