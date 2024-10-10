#!/usr/bin/env Rscript

writeLines("\n### Running _1_get_ctss_from_bw.r ###")

writeLines("\n# Importing R libraries..")
suppressPackageStartupMessages({
  library(argparse)
  library(CAGEfightR)
  library(GenomicRanges)
  library(PRIMEloci)
})



### ARGPARSE
parser <- ArgumentParser()

# Input
parser$add_argument("-i", "--input_dir", default = "./",
                    help = "CAGE BigWig directory")
parser$add_argument("-m", "--design_matrix", default = "design.matrix.tsv",
                    help = "Design matrix of all CAGE data in .tsv format, example can be found in the 'example/resources' directory") # nolint: line_length_linter.

# Output
parser$add_argument("-o", "--output_dir", default = "./",
                    help = "Output directory")
parser$add_argument("-n", "--outfile", default = "ctss_rse.rds",
                    help = "Output file name for ctss_rse object")

# Parameters
parser$add_argument("-k", "--keep_standard_chr", action = "store_true",
                    help = "Keep only standard chromosomes")

args <- parser$parse_args()

# Setting up variables
input_dir <- args$input_dir
design_matrix_file <- args$design_matrix

output_dir <- args$output_dir
outfile_ctss_rse <- args$outfile



# Read in CAGE BigWig, based on design matrix provided,
# no filtering/subsetBySupport steps are performed
writeLines("\n# Reading in data..")
design_matrix <- read.table(design_matrix_file,
                            header = TRUE,
                            sep = "\t",
                            row.names = "Name")
writeLines("\n# Getting CTSSs from BigWig files..")
ctss_rse <- get_ctss_from_bw(input_dir, design_matrix)

# Keep only standard chromosomes if specified
if (args$keep_standard_chr) {
  writeLines("\n# Keeping only standard chromosomes..")
  ctss_rse <- GenomeInfoDb::keepStandardChromosomes(ctss_rse,
                                                    pruning.mode = "coarse")
}

# Save
writeLines("\n# Saving ctss object..\n")
saveRDS(ctss_rse, file.path(output_dir, outfile_ctss_rse))