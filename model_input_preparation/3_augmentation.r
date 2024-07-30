#!/usr/bin/env Rscript

writeLines("\n #### 3 DATA AUGMENTATION ####")

writeLines("\nImporting R libraries..")
suppressPackageStartupMessages({
  library(tools)
  library(argparse)
  library(tidyverse)
  library(dplyr)
  library(CAGEfightR)
  library(AnnotationDbi)
  library(GenomicFeatures)
  library(GenomicRanges)
  library(PRIMEloci)
})

set.seed(214)

### ARGPARSE
parser <- ArgumentParser()
parser$add_argument("-r", "--ocr_rdata", default = "./",
                    help = "OCR and OCR-like RData file")
parser$add_argument("-n", "--name", default = "name", help = "Name")
parser$add_argument("-o", "--output_dir", default = "./",
                    help = "output directory")
parser$add_argument("-d", "--distance", default = 200, type = "integer",
                    help = "distance from the DHS peak (bps)")

args <- parser$parse_args()

# general
dist <- args$distance

# input and output
ocr_rdata_file <- args$ocr_rdata
output_dir <- args$output_dir
name <- args$name

# 1 load Rdata
load(ocr_rdata_file)

# 2 augment the positive OCR based on
# normal distribution prob within +- 30 bps from center
# they were check to make sure that
# the augment one not exact overlap with original
writeLines("\naugmenting the training- positive- ocr based on normal distribution prob within +- 30 bps from center..") # nolint: line_length_linter.
ocr_train_aug_gr <- augmentation_norm(ocr_train_gr)
print(paste0("train positive: ", length(ocr_train_gr)))
print(paste0("augmented train positive: ", length(ocr_train_aug_gr)))

# 3 augment the negative OCR based on
# random distribution prob within +- 30 bps from center
# they were check to make sure tha
# the augment one not exact overlap with original
writeLines("\naugmenting the trianing- negative- ocr based on random distribution prob within +- 30 bps from center..") # nolint: line_length_linter.
ocrlike_neg_train_aug_gr <- augmentation_unif(ocrlike_neg_train_gr)
print(paste0("train negative: ", length(ocrlike_neg_train_gr)))
print(paste0("augmented train negative (before filtering): ",
             length(ocrlike_neg_train_aug_gr)))

# 4 subset the training negative to training positive (orig + aug)
# we allow the negative to be "50% or less" overlap with the positive
ocrlike_neg_train_aug_rm <- IRanges::subsetByOverlaps(ocrlike_neg_train_aug_gr,
                                                      c(ocr_train_aug_gr,
                                                        ocr_train_gr),
                                                      minoverlap = dist)
print(paste0("remove augmented train negative that overlap more than 50% overlap with pos): ", # nolint: line_length_linter.
             length(ocrlike_neg_train_aug_rm)))
ocrlike_neg_train_aug_gr <- IRanges::subsetByOverlaps(ocrlike_neg_train_aug_gr,
                                                      c(ocr_train_aug_gr,
                                                        ocr_train_gr),
                                                      minoverlap = dist,
                                                      invert = TRUE)
print(paste0("augmented train negative (remove more than 50% overlap with pos): ", # nolint: line_length_linter.
             length(ocrlike_neg_train_aug_gr)))

# 5 save
writeLines("\nSaving..")
save(list = c("ocr_train_gr",
              "ocr_train_aug_gr",
              "ocr_test_gr",
              "ocrlike_neg_train_gr",
              "ocrlike_neg_train_aug_gr",
              "ocrlike_neg_test_gr"),
     file = paste0("3_", name, "_ocr_ocrlikeneg_augmentation.RData"))
