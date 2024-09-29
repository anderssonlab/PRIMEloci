#!/usr/bin/env Rscript

writeLines("\n### Running get_tc_profiles.r ###")

writeLines("\n# Importing R libraries..")
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(parallel)
  library(argparse)
})

# Create argument parser
parser <- ArgumentParser(description = "Sliding window on genomic ranges")
parser$add_argument("-i", "--input", required = TRUE,
                    help = "Path to the input .rds file (tc_grl)")
parser$add_argument("-o", "--output_dir", default = "./",
                    help = "Output directory")
parser$add_argument("-t", "--outfile_sld_tc_grl", default = "sld_tc_grl.RDS",
                    help = "Output file name for sld_tc_grl object")
parser$add_argument("-s", "--slide_by", type = "integer", default = 20,
                    help = "Slide by parameter")
parser$add_argument("-e", "--expand_by", type = "integer", default = 200,
                    help = "Expand by parameter")

# Parse arguments
args <- parser$parse_args()

# Load the input data
tc_grl <- readRDS(args$input)

# Output
output_dir <- args$output_dir
outfile_sld_tc_grl <- args$outfile_sld_tc_grl

# Log the current time
current_datetime <- Sys.time()
formatted_datetime <- format(current_datetime, "%Y-%m-%d %H:%M:%S")
cat("Start time:", formatted_datetime, "\n")

# Sliding window operation
genomewide_sliding_window_grl <- lapply(tc_grl, function(gr) {
  genomewide_sliding_window(gr,
                            slide_by = args$slide_by,
                            expand_by = args$expand_by,
                            use_max_cores = TRUE)
})

# Convert to GRangesList
genomewide_sliding_window_grl <- GenomicRanges::GRangesList(genomewide_sliding_window_grl) # nolint: line_length_linter.

# Log the current time after processing
current_datetime <- Sys.time()
formatted_datetime <- format(current_datetime, "%Y-%m-%d %H:%M:%S")
cat("End time:", formatted_datetime, "\n")

# Save the result to the specified output file
writeLines("\n# Saving tc objects..\n")
saveRDS(genomewide_sliding_window_grl,
        file.path(output_dir, outfile_sld_tc_grl))