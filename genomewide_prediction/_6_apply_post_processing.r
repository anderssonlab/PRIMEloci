#!/usr/bin/env Rscript

writeLines("\n### Running _filter_core_overlapping.r ###")

writeLines("\n# Importing R libraries..")
suppressPackageStartupMessages({
  library(argparse)
  library(GenomicRanges)
  library(parallel)
  library(tools)
  library(stringr)
  library(PRIMEloci)
})

# Create the argparse object
parser <- ArgumentParser(description = "Process GRanges based on a BED file and score threshold.") # nolint: line_length_linter.
parser$add_argument("-i", "--input_file", required = TRUE,
                    help = "Path to the input BED file.")
parser$add_argument("-t", "--score_threshold", type = "double", default = 0.7,
                    help = "Score threshold for filtering GRanges.")
parser$add_argument("-c", "--use_max_cores", action = "store_true",
                    default = TRUE,
                    help = "Flag to use the maximum number of available cores.")
parser$add_argument("-o", "--output_dir", default = NULL,
                    help = "Directory to save the output files.")

# Parse arguments
args <- parser$parse_args()

# Define parameters based on parsed arguments
score_threshold <- args$score_threshold
bed_file <- args$input_file
use_max_cores <- args$use_max_cores
output_dir <- args$output_dir

# Set the number of cores
num_cores <- if (use_max_cores) { parallel::detectCores() } else { 1 }

# Load and prepare data
bed <- load_bed_file(bed_file)
gr <- create_granges_from_bed(bed)

# Filter GRanges by score threshold
filtered_gr <- gr[gr$score >= score_threshold]

# Get a list of unique chromosomes
chr_list <- unique(as.character(seqnames(filtered_gr)))

# Run in parallel across chromosomes using the specified number of cores
collapsed_gr_list <- mclapply(chr_list,
                              process_by_chr,
                              filtered_gr = filtered_gr,
                              mc.cores = num_cores)

# Combine all results
collapsed_gr <- do.call(c, collapsed_gr_list)
collapsed_gr <- sort(collapsed_gr)

# If output_dir is specified, save the GRanges object to the directory
if (!is.null(output_dir) && output_dir != FALSE) {
  input_basename <- tools::file_path_sans_ext(basename(bed_file))
  input_basename <- input_basename %>%
    stringr::str_replace_all("all", as.character(score_threshold)) %>%  # Replace "all" with threshold # nolint: line_length_linter.
    stringr::str_replace_all("[^[:alnum:]]", "_")                # Replace non-alphanumeric characters with "_" # nolint: line_length_linter.

  save_granges_to_bed_2(collapsed_gr, output_dir, input_basename)
}

writeLines("Done!")