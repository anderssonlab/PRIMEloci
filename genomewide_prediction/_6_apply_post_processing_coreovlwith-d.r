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

# Create the argparse object
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



# Parse arguments
args <- parser$parse_args()

# Define parameters based on parsed arguments
score_threshold <- args$score_threshold
score_diff <- args$score_diff
bed_file <- args$input_file
use_max_cores <- args$use_max_cores
output_dir <- args$output_dir

core_width <- 151

# Set the number of cores
num_cores <- parallel::detectCores()
if (!use_max_cores) {
  num_cores <- 1
}

# Load and prepare data
bed <- load_bed_file(bed_file)
gr <- create_granges_from_bed(bed)

# Filter GRanges by score threshold
filtered_gr <- gr[gr$score >= score_threshold]
# Get a list of unique chromosomes
chr_list <- unique(as.character(seqnames(filtered_gr)))



collapsed_gr_list <- mclapply(chr_list, function(chr) {

  tryCatch({

    chr_gr <- filtered_gr[GenomicRanges::seqnames(filtered_gr) == chr]

    core_gr <- GenomicRanges::resize(chr_gr, width = core_width, fix = "center")
    core_gr$thick <- start(core_gr) + floor(core_width / 2)

    final_merged_gr <- selective_merge_cores(core_gr, score_diff)

    return(final_merged_gr)
  }, error = function(e) {
    message(paste("Error processing chromosome", chr, ":", e$message))
    return(NULL)
  })
}, mc.cores = num_cores)

collapsed_gr <- do.call(c, collapsed_gr_list)
collapsed_gr <- sort(collapsed_gr)

if (!is.null(output_dir) && output_dir != FALSE) {
  input_basename <- tools::file_path_sans_ext(basename(bed_file))
  input_basename <- input_basename %>%
    stringr::str_replace_all("all", as.character(score_threshold)) %>%
    stringr::str_replace_all("[^[:alnum:]]", "_")
  write_granges_to_bed_coreovlwithd(collapsed_gr,
                                    output_dir,
                                    input_basename,
                                    score_diff)
}

writeLines("Done!+++")