writeLines("\n### Running _4_get_profile.r ###")

writeLines("\nImporting R libraries..")
suppressPackageStartupMessages({
  library(argparse)
  library(CAGEfightR)
  library(parallel)
  library(GenomicRanges)
  library(PRIMEloci)
})


# Define argument parser
parser <- ArgumentParser()

# Add arguments
# Input
parser$add_argument("-c", "--infile_ctss_rse", default = "./ctss_rse.RDS",
                    help = "FULL PATH to input file name for ctss_rse rds object") # nolint: line_length_linter.
parser$add_argument("-t", "--infile_tc_grl", default = "./tc_grl.RDS",
                    help = "FULL PATH to input file name for tc_grl rds object")

# Output
parser$add_argument("-o", "--output_dir", default = "./",
                    help = "Profile output directory")
parser$add_argument("-n", "--name", default = "name",
                    help = "Name of main profile directory")
parser$add_argument("-s", "--save_count_profiles", action = "store_true",
                    default = FALSE,
                    help = "Flag to save count profile. Default is FALSE.")
parser$add_argument("-f", "--file_format", type = "character",
                    default = "parquet",
                    choices = c("parquet", "csv"),
                    help = "File format for output files. Choose between 'parquet' or 'csv'. Default is 'parquet'.") # nolint: line_length_linter.

# Parameters
parser$add_argument("-e", "--ext_dis", default = 200,
                    help = "Extension distance")

# Parse arguments
args <- parser$parse_args()

# Setting up variables
infile_ctss_rse <- args$infile_ctss_rse
infile_tc_grl <- args$infile_tc_grl

outdir_dir <- args$output_dir
outdir_dir_name <- args$name
file_format <- args$file_format
save_count_profiles <- args$save_count_profiles

ext_dis <- as.numeric(args$ext_dis)


# Read in .RDS files
writeLines("\nReading in data..")
ctss_rse <- readRDS(infile_ctss_rse)
tc_grl <- readRDS(infile_tc_grl)


# Create profiles
writeLines(paste0("\nCreating profile.."))
report_time_execution(PRIMEloci_profile(ctss_rse,
                                        tc_grl,
                                        outdir_dir,
                                        outdir_dir_name,
                                        ext_dis,
                                        save_count_profiles = save_count_profiles, # nolint: line_length_linter.
                                        file_type = file_format))

writeLines("\n### Finished _4_get_profile.r ###\n")
