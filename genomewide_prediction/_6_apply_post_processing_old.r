#!/usr/bin/env Rscript

writeLines("\n### Running _6_apply_post_processing.r ###")

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

library(GenomicRanges)
library(S4Vectors)
library(IRanges)
library(dplyr)
ovlcore_reduced_by_chr_2 <- function(chr, filtered_gr, score_diff) {
  tryCatch({
    # Subset GRanges by chromosome
    chr_gr <- filtered_gr[GenomicRanges::seqnames(filtered_gr) == chr]
    core_gr <- GenomicRanges::resize(chr_gr, width = 151, fix = "center")

    #overlaps <- GenomicRanges::findOverlaps(core_gr)
    #overlaps <- overlaps[overlaps$score >= (max(overlaps$score) - score_diff)]
    #collapsed_ranges <- GenomicRanges::reduce(core_gr[unique(S4Vectors::queryHits(overlaps))]) # nolint: line_length_linter.

    # Initialize metadata columns
    thick_vals <- numeric(length(collapsed_ranges))
    max_scores <- numeric(length(collapsed_ranges))
    all_scores <- character(length(collapsed_ranges))

    # Loop through each collapsed range to extract metadata
    for (i in seq_along(collapsed_ranges)) {
      metadata <- get_metadata(collapsed_ranges[i], filtered_gr)
      thick_vals[i] <- metadata$thick
      max_scores[i] <- metadata$max_score
      all_scores[i] <- metadata$all_scores
    }

    # Assign metadata to collapsed ranges
    collapsed_ranges$thick <- thick_vals
    collapsed_ranges$max_score <- max_scores
    collapsed_ranges$all_scores <- all_scores

    return(collapsed_ranges)
  }, error = function(e) {
    message(paste("Error processing chromosome", chr, ":", e$message))
    return(NULL)  # Return NULL on error
  })
}


# Load and prepare data
bed <- load_bed_file(bed_file)
gr <- create_granges_from_bed(bed)

# Filter GRanges by score threshold
filtered_gr <- gr[gr$score >= score_threshold]
# Get a list of unique chromosomes
chr_list <- unique(as.character(seqnames(filtered_gr)))

# Run in parallel across chromosomes using the specified number of cores
collapsed_gr_list <- mclapply(chr_list,
                              ovlcore_reduced_by_chr_2,
                              filtered_gr = filtered_gr,
                              score_diff = 0.15,
                              mc.cores = num_cores)
print(collapsed_gr_list)

# Combine all results
collapsed_gr <- do.call(c, collapsed_gr_list)
collapsed_gr <- sort(collapsed_gr)

# If output_dir is specified, save the GRanges object to the directory
if (!is.null(output_dir) && output_dir != FALSE) {
  input_basename <- tools::file_path_sans_ext(basename(bed_file))
  input_basename <- input_basename %>%
    stringr::str_replace_all("all", as.character(score_threshold)) %>%  # Replace "all" with threshold # nolint: line_length_linter.
    stringr::str_replace_all("[^[:alnum:]]", "_")                # Replace non-alphanumeric characters with "_" # nolint: line_length_linter.

  save_ovlcorereduced_to_bed(collapsed_gr, output_dir, input_basename)
}

writeLines("Done!")
