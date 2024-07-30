# Load necessary libraries
suppressPackageStartupMessages({

  library(dotenv)
  library(this.path)

  library(tidyverse)

  library(IRanges)
  library(GenomeInfoDb)
  library(GenomicRanges)
  library(GenomicFeatures)

  library(PRIMEloci)

})

# Working directory and environment variables
work_dir <- this.path::this.dir()
setwd(work_dir)
dotenv::load_dot_env("epcrispr_benchmark.env")


# 0. setting

# Retrieve environment variables
ext_dis <- as.numeric(Sys.getenv("EXT_DIS"))

# main directoriesi
dir_resources <- Sys.getenv("DIR_RESOURCES")
dir_results <- Sys.getenv("DIR_RESULTS")

# epc profiles directory's names
outdir_dir <- Sys.getenv("OUTPUT_DIR_EPCPROFILES")
outdir_dir_name <- Sys.getenv("OUTPUT_DIR_NAME_EPCPROFILES")

# set this according to the data
outdir_main_name <- c("metadata",
                      "profiles",
                      "profiles_subtnorm",
                      "predictions")
outdir_subdir_name <- c("pos", "neg")
outfile_pos_gr <- Sys.getenv("OUTFILE_POS_GR")
outfile_neg_notsig_gr <- Sys.getenv("OUTFILE_NEG_NOTSIG_GR")

# read in ranges data
pos_gr <- readRDS(file.path(dir_results, outfile_pos_gr))
neg_gr <- readRDS(file.path(dir_results, outfile_neg_notsig_gr))

# read in CTSSs data
infile_ctss_rse <- Sys.getenv("INFILE_CTSS_RSE")
ctss_rse <- readRDS(file.path(dir_resources, infile_ctss_rse))

# 1 create the output dir
prep_profile_dir(output_dir = outdir_dir,
                 output_dir_name = outdir_dir_name,
                 output_main_name = outdir_main_name,
                 output_subdir_name = outdir_subdir_name)

# 2. create profile
# modify region_gr and outdir_subdir_name
report_time_execution(wrapup_make_profiles(ctss_rse, pos_gr,
                                           dir_results, outdir_dir_name,
                                           "pos", ext_dis,
                                           addtn_to_filename = "_EPC_pos",
                                           save_count_profiles = TRUE))
report_time_execution(wrapup_make_profiles(ctss_rse, neg_gr,
                                           dir_results, outdir_dir_name,
                                           "neg", ext_dis,
                                           addtn_to_filename = "_EPC_neg",
                                           save_count_profiles = TRUE))
