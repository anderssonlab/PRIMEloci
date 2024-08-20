#!/usr/bin/env Rscript

writeLines("\n### Running _filter_bed_to_reduce_gr.r ###")

writeLines("\n# Importing R libraries..")
suppressPackageStartupMessages({
  library(assertthat)
  library(data.table)
  library(argparse)
  library(foreach)
  library(doParallel)
  library(GenomicRanges)
  library(rtracklayer)
  library(PRIMEloci)
})

# Define the command line arguments
parser <- ArgumentParser()
parser$add_argument("-i", "--input_bed", type = "character",
                    help = "Input BED file")
parser$add_argument("-o", "--output_dir", type = "character",
                    default = NULL,
                    help = "Output directory, default will be the same as the input file")
args <- parser$parse_args()



# Execute the main function
wrapup_filter_bed_to_reduce(args$input_bed, args$output_dir)

writeLines("Done!")
