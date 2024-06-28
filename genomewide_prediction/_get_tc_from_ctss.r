#!/usr/bin/env Rscript

writeLines("\n### Running get_tc_profiles.r ###")

writeLines("\n# Importing R libraries..")
suppressPackageStartupMessages({
  library(argparse)
  library(CAGEfightR)
  library(GenomicRanges)
})

source("../R/tc.r")

### ARGPARSE
parser <- ArgumentParser()

# Input
parser$add_argument("-c", "--infile_ctss_rse", default = "./ctss_rse.RDS",
                    help = "FULL PATH to input file name for ctss_rse rds object") # nolint: line_length_linter.

# Output
parser$add_argument("-o", "--output_dir", default = "./",
                    help = "Output directory")
parser$add_argument("-t", "--outfile_tc_grl", default = "tc_grl.RDS",
                    help = "Output file name for tc_grl object")

# Parameters
parser$add_argument("-e", "--ext_dis", default = 200,
                    help = "Extension distance")

args <- parser$parse_args()

# Setting up variables
infile_ctss_rse <- args$infile_ctss_rse
ctss_rse <- readRDS(infile_ctss_rse)

output_dir <- args$output_dir
outfile_tc_grl <- args$outfile_tc_grl

ext_dis <- as.integer(args$ext_dis)

writeLines("\n# Gettign tag clusters and extending the start/end using thick..\n") # nolint: line_length_linter.
tc_grl <- get_tcs_and_extend_fromthick(ctss_rse, ext_dis = ext_dis)

# Save
writeLines("\n# Saving tc objects..\n")
saveRDS(tc_grl, file.path(output_dir, outfile_tc_grl))