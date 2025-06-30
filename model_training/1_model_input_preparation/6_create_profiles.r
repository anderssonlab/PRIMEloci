#!/usr/bin/env Rscript

writeLines("\n #### 6 CREATE PROFILES ####")

writeLines("\nImporting R libraries..")
suppressPackageStartupMessages({
  library(argparse)
  library(tidyverse)
  library(dplyr)
  library(CAGEfightR)
  library(AnnotationDbi)
  library(GenomicFeatures)

  library(rtracklayer)
  library(tibble)
  library(readr)
  library(forcats)
  library(tidyr)
  library(GenomicRanges)
  library(reshape2)

  library(PRIMEloci)
})

### ARGPARSE
parser <- ArgumentParser()
parser$add_argument("-r", "--ocr_rdata", default = ".",
                    help = "OCR and OCR-like RData file ")
parser$add_argument("-c", "--cage_tss_rdata", default = ".",
                    help = "CAGE TSS rData file")
parser$add_argument("-o", "--output_dir", default = ".",
                    help = "output directory")
parser$add_argument("-O", "--output_dir_name", default = "profiles_result",
                    help = "output dir name")
parser$add_argument("-d", "--distance", default = 200, type = "integer",
                    help = "distance from the dhs peak (bps)")
parser$add_argument("-s", "--save_count_profiles", action = "store_true",
                    default = FALSE,
                    help = "Flag to save count profile. Default is FALSE.")

args <- parser$parse_args()

# input and output
ocr_rdata_file <- args$ocr_rdata
tss_rdata_file <- args$cage_tss_rdata
outdir_dir <- args$output_dir
outdir_dir_name <- args$output_dir_name

# Fixed
ext_dis <- as.numeric(args$distance)
save_count_profiles <- args$save_count_profiles

outdir_main_name <- c("metadata",
                      "profiles",
                      "profiles_subtnorm",
                      "predictions")
outdir_subdir_name <- c("pos", "neg")

# 1 load Rdata
load(ocr_rdata_file)
load(tss_rdata_file)

# 2 create the output dir
prep_profile_dir(output_dir = outdir_dir,
                 output_dir_name = outdir_dir_name,
                 output_main_name = outdir_main_name,
                 output_subdir_name = outdir_subdir_name)

# 3 create profile and its metadata

# Define a function to create profiles for a single subdir_name
report_time_execution(wrapup_make_profiles(tr_orig_ctss, tr_vl_ls$pos_tr_grl,
                                           outdir_dir, outdir_dir_name,
                                           "pos", ext_dis,
                                           addtn_to_filename = "_tr",
                                           save_count_profiles = TRUE,
                                           file_type = "csv"))
report_time_execution(wrapup_make_profiles(tr_orig_ctss, tr_vl_ls$pos_vl_grl,
                                           outdir_dir, outdir_dir_name,
                                           "pos", ext_dis,
                                           addtn_to_filename = "_vl",
                                           save_count_profiles = TRUE,
                                           file_type = "csv"))

writeLines("DONE POSITIVE..")

report_time_execution(wrapup_make_profiles(tr_orig_ctss, tr_vl_ls$neg_tr_grl,
                                           outdir_dir, outdir_dir_name,
                                           "neg", ext_dis,
                                           addtn_to_filename = "_tr",
                                           save_count_profiles = TRUE,
                                           file_type = "csv"))
report_time_execution(wrapup_make_profiles(tr_orig_ctss, tr_vl_ls$neg_vl_grl,
                                           outdir_dir, outdir_dir_name,
                                           "neg", ext_dis,
                                           addtn_to_filename = "_vl",
                                           save_count_profiles = TRUE,
                                           file_type = "csv"))

writeLines("DONE NEGATIVE..")

report_time_execution(wrapup_make_profiles(te_orig_ctss, te_ls$pos_te_grl,
                                           outdir_dir, outdir_dir_name,
                                           "pos", ext_dis,
                                           addtn_to_filename = "_te",
                                           save_count_profiles = TRUE,
                                           file_type = "csv"))
report_time_execution(wrapup_make_profiles(te_orig_ctss, te_ls$neg_te_grl,
                                           outdir_dir, outdir_dir_name,
                                           "neg", ext_dis,
                                           addtn_to_filename = "_te",
                                           save_count_profiles = TRUE,
                                           file_type = "csv"))

report_time_execution(wrapup_make_profiles(te_orig_ctss, te_aug_ls$pos_te_grl,
                                           outdir_dir, outdir_dir_name,
                                           "pos", ext_dis,
                                           addtn_to_filename = "_te_aug",
                                           save_count_profiles = TRUE,
                                           file_type = "csv"))
report_time_execution(wrapup_make_profiles(te_orig_ctss, te_aug_ls$neg_te_grl,
                                           outdir_dir, outdir_dir_name,
                                           "neg", ext_dis,
                                           addtn_to_filename = "_te_aug",
                                           save_count_profiles = TRUE,
                                           file_type = "csv"))

writeLines("DONE TESTING..")

report_time_execution(wrapup_make_profiles(tr_orig_ctss,
                                           tr_vl_aug_ls$pos_tr_grl,
                                           outdir_dir, outdir_dir_name,
                                           "pos", ext_dis,
                                           addtn_to_filename = "_tr_aug",
                                           save_count_profiles = TRUE,
                                           file_type = "csv"))
report_time_execution(wrapup_make_profiles(tr_orig_ctss,
                                           tr_vl_aug_ls$pos_vl_grl,
                                           outdir_dir, outdir_dir_name,
                                           "pos", ext_dis,
                                           addtn_to_filename = "_vl_aug",
                                           save_count_profiles = TRUE,
                                           file_type = "csv"))
report_time_execution(wrapup_make_profiles(tr_orig_ctss,
                                           tr_vl_aug_ls$neg_tr_grl,
                                           outdir_dir, outdir_dir_name,
                                           "neg", ext_dis,
                                           addtn_to_filename = "_tr_aug",
                                           save_count_profiles = TRUE,
                                           file_type = "csv"))
report_time_execution(wrapup_make_profiles(tr_orig_ctss,
                                           tr_vl_aug_ls$neg_vl_grl,
                                           outdir_dir, outdir_dir_name,
                                           "neg", ext_dis,
                                           addtn_to_filename = "_vl_aug",
                                           save_count_profiles = TRUE,
                                           file_type = "csv"))
writeLines("DONE AUGMENTATION..")