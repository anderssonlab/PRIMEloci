#!/usr/bin/env Rscript

writeLines("\n #### 5 VALIDATION AND SELECTED NEGATIVE ####")

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

set.seed(214)

### ARGPARSE
parser <- ArgumentParser()
parser$add_argument("-r", "--ocr_rdata", default = "./",
                    help = "OCR and OCR-like RData file ")
parser$add_argument("-c", "--cage_tss_rdata", default = ".",
                    help = "CAGE TSS rData file")
parser$add_argument("-n", "--name", default = "name",
                    help = "Name")
parser$add_argument("-v", "--validation_ratio", default = 0.2,
                    help = "validation_ratio")

args <- parser$parse_args()

# input and output
ocr_rdata_file <- args$ocr_rdata
tss_rdata_file <- args$cage_tss_rdata
name <- args$name
validation_ratio <- args$validation_ratio



# 1 load Rdata
writeLines("\nLoad data..")
load(ocr_rdata_file)
load(tss_rdata_file)


# 2 select OCR for each library with condition
#In summary, the function filters OCRs based on two main conditions:
# 1. OCRs with more than one-base overlap are always selected.
# 2. Singleton OCRs are selected only if their score is
# greater than or equal to min_score_for_onebase_ocr (default is 4).
writeLines("\nSelect OCR for each library..") # nolint: line_length_linter.
writeLines("\n..")

tr_pos_sltocr_grl <- select_ocr_on_condition(tr_orig_ctss, ocr_train_gr)
tr_pos_aug_sltocr_grl <- select_ocr_on_condition(tr_orig_ctss, ocr_train_aug_gr)

te_pos_sltocr_grl <- select_ocr_on_condition(te_orig_ctss, ocr_test_gr)
te_pos_aug_sltocr_grl <- select_ocr_on_condition(te_orig_ctss, ocr_test_aug_gr)

# tr_neg_sltocr_grl <- select_ocr_on_condition(tr_orig_ctss, ocrlike_neg_train_gr) # nolint: line_length_linter.

tr_neg_aug_sltocr_grl <- select_ocr_on_condition(tr_orig_ctss, ocrlike_neg_train_aug_gr) # nolint: line_length_linter.
tr_neg_npn_sltocr_grl <- select_ocr_on_condition(tr_orig_ctss, nearpos_neg_train_aug_gr) # nolint: line_length_linter.

te_neg_aug_sltocr_grl <- select_ocr_on_condition(te_orig_ctss, ocrlike_neg_test_aug_gr) # nolint: line_length_linter.
te_neg_npn_sltocr_grl <- select_ocr_on_condition(te_orig_ctss, nearpos_neg_test_aug_gr) # nolint: line_length_linter.

# te_neg_sltocr_grl <- select_ocr_on_condition(te_orig_ctss, ocrlike_neg_test_gr) # nolint: line_length_linter.



# 3 split data
writeLines("\nSplit data..")

# TRAINING POSITIVE + TESTING POSITIVE
tr_pos_grl <- subset_annotation_positive(tr_pos_sltocr_grl)
tr_pos_aug_grl <- subset_annotation_positive(tr_pos_aug_sltocr_grl)
te_pos_grl <- subset_annotation_positive(te_pos_sltocr_grl)
te_pos_aug_grl <- subset_annotation_positive(te_pos_aug_sltocr_grl)

# TRAINING NEGATIVE + TESTING NEGATIVE
# subsampling negative OCR position to match number of positive,
# but keep annotation ratio of negative
# nearpos negatives were match to positive and
# aug negatives were match to aug positive
tr_neg_npn_grl <- select_negative_to_match_positive(tr_neg_npn_sltocr_grl, tr_pos_grl) # nolint: line_length_linter.
tr_neg_aug_grl <- select_negative_to_match_positive(tr_neg_aug_sltocr_grl, tr_pos_aug_grl) # nolint: line_length_linter.
te_neg_npn_grl <- select_negative_to_match_positive(te_neg_npn_sltocr_grl, te_pos_grl) # nolint: line_length_linter.
te_neg_aug_grl <- select_negative_to_match_positive(te_neg_aug_sltocr_grl, te_pos_aug_grl) # nolint: line_length_linter.




# TRAINING: split training (80%) and validation (20%)
# based on their natural ratio of annotation
tr_vl_ls <- wrapup_split_train_valid(tr_pos_grl, tr_neg_npn_grl)
tr_vl_aug_ls <- wrapup_split_train_valid(tr_pos_aug_grl, tr_neg_aug_grl)



# TESTING: keep equal number
wrapup_test <- function(pos_grl, neg_grl) {
  pos_te_grl <- GenomicRanges::GRangesList()
  neg_te_grl <- GenomicRanges::GRangesList()

  for (i in seq_along(pos_grl)) {
    final_num_te <- min(length(pos_grl[[i]]), length(neg_grl[[i]]))

    pos_te_grl[[i]] <- keep_equal_num(pos_grl[[i]], final_num_te)
    neg_te_grl[[i]] <- keep_equal_num(neg_grl[[i]], final_num_te)
  }

  # make thick column to IRanges
  pos_te_grl <- convert_to_thick_iranges(pos_te_grl)
  neg_te_grl <- convert_to_thick_iranges(neg_te_grl)

  te_ls <- list(pos_te_grl = pos_te_grl,
                neg_te_grl = neg_te_grl)

  return(te_ls)

}

te_ls <- wrapup_test(te_pos_grl, te_neg_npn_grl)
te_aug_ls <- wrapup_test(te_pos_aug_grl, te_neg_aug_grl)


# 4 save
writeLines("\nSaving..")
save(list = c("tr_vl_ls", "tr_vl_aug_ls", "te_ls", "te_aug_ls"),
     file = paste0("5_", name, "_tr_vl_te.RData"))