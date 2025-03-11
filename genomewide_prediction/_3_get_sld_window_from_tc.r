#!/usr/bin/env Rscript

writeLines("\n### Running _3_get_sld_window_from_tc.r ###")

writeLines("\n# Importing R libraries..")
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(parallel)
  library(argparse)
  library(PRIMEloci)
})

# Create argument parser
parser <- ArgumentParser(description = "Sliding window on genomic ranges")
parser$add_argument("-i", "--infile", required = TRUE,
                    help = "Path to the input .rds file (tc_grl)")

parser$add_argument("-o", "--output_dir", default = "./",
                    help = "Output directory")
parser$add_argument("-n", "--outfile", default = "sld_tc_grl.rds",
                    help = "Output file name for sld_tc_grl object")

parser$add_argument("-s", "--sld_by", type = "integer", default = 20,
                    help = "Slide by parameter")
parser$add_argument("-e", "--ext_dis", default = 200,
                    help = "Extension distance")

# Parse arguments
args <- parser$parse_args()

# Load the input data
tc_grl <- readRDS(args$infile)

# Output
output_dir <- args$output_dir
outfile_sld_tc_grl <- args$outfile

# Log the current time
current_datetime <- Sys.time()
formatted_datetime <- format(current_datetime, "%Y-%m-%d %H:%M:%S")
cat("Start time:", formatted_datetime, "\n")

# Sliding window operation
tc_sliding_window_grl <- lapply(seq_along(tc_grl), function(i) {
  print("Sliding window operation..")

  # Print the name of the GRanges object
  gr_name <- names(tc_grl)[i]  # Extract name from the list level
  print(paste("Processing:", gr_name))

  # Get the actual GRanges object
  gr <- tc_grl[[i]]

  # Run sliding window operation
  tc_sliding_window(gr,
                    slide_by = as.numeric(args$sld_by),
                    expand_by = as.numeric(args$ext_dis),
                    num_cores = NULL)
})

# Convert to GRangesList
tc_sliding_window_grl <- GenomicRanges::GRangesList(tc_sliding_window_grl) # nolint: line_length_linter.

# Log the current time after processing
current_datetime <- Sys.time()
formatted_datetime <- format(current_datetime, "%Y-%m-%d %H:%M:%S")
cat("End time:", formatted_datetime, "\n")

# Save the result to the specified output file
writeLines("\n# Saving tc objects..\n")
saveRDS(tc_sliding_window_grl, file.path(output_dir, outfile_sld_tc_grl))