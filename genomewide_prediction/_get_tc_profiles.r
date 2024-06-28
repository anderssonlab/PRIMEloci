writeLines("\n### Running get_tc_profiles.r ###")

writeLines("\nImporting R libraries..")
suppressPackageStartupMessages({
  library(argparse)
  library(CAGEfightR)
  library(GenomicRanges)
})

# Source additional script
source("../R/general.r")
source("../R/manipulate_gr.r")
source("../R/directories.r")
source("../R/profiles.r")
source("../R/heatmap.r")

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
parser$add_argument("-r", "--output_subdir_name",
                    type = "character", default = "tcs",
                    help = "Comma-separated list of subdirectory names")
parser$add_argument("-s", "--save_count_profiles", action = "store_true",
                    default = FALSE,
                    help = "Flag to save count profile. Default is FALSE.")

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
save_count_profiles <- args$save_count_profiles

ext_dis <- as.numeric(args$ext_dis)

outdir_main_name <- c("metadata",
                      "profiles",
                      "profiles_subtnorm",
                      "predictions")
outdir_subdir_name <- unlist(strsplit(args$output_subdir_name, ","))


# Read in .RDS files
writeLines("\nReading in data..")
ctss_rse <- readRDS(infile_ctss_rse)
tc_grl <- readRDS(infile_tc_grl)


# Create output directory
prep_profile_dir(output_dir = outdir_dir,
                 output_dir_name = outdir_dir_name,
                 output_main_name = outdir_main_name,
                 output_subdir_name = outdir_subdir_name)


# Define a function to create profiles for a single subdir_name
create_profiles <- function(subdir_name) {
  writeLines(paste0("\nCreating profiles for ", subdir_name, ".."))
  report_time_execution(wrapup_make_profiles(ctss_rse,
                                             tc_grl,
                                             outdir_dir,
                                             outdir_dir_name,
                                             subdir_name,
                                             ext_dis,
                                             save_count_profiles = save_count_profiles)) # nolint: line_length_linter.
}

lapply(outdir_subdir_name, create_profiles)

writeLines("\n### Finished get_tc_profiles.r ###\n")