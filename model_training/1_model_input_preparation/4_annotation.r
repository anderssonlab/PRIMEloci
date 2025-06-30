#!/usr/bin/env Rscript

writeLines("\n #### 4 ANNOTATION ####")

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

# need to be changed if changing the genome
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

### ARGPARSE
parser <- ArgumentParser()
parser$add_argument("-r", "--ocr_rdata", default = ".",
                    help = "OCR and OCR-like RData file ")
parser$add_argument("-c", "--cage_tss_rdata", default = "./",
                    help = "CAGE TSS rData file")
parser$add_argument("-n", "--name", default = "name",
                    help = "Name")
parser$add_argument("-o", "--output_dir", default = "./",
                    help = "output directory")
parser$add_argument("-u", "--ucsc", default = "TRUE",
                    help = "USE UCSC hg38 for txdb")
parser$add_argument("-a", "--annotation_file",
                    default = "/projects/ralab/data/projects/nucleiCAGEproject/0.External_resources/gencode.v31.annotation.gtf", # nolint: line_length_linter.
                    help = "gtf annotation file from GENCODE")
parser$add_argument("-g", "--ucscgenome", default = "hg38",
                    help = "UCSC genome")
parser$add_argument("-s", "--species", default = "Homo sapiens",
                    help = "species")

args <- parser$parse_args()

# input and output
ocr_rdata_file <- args$ocr_rdata
tss_rdata_file <- args$cage_tss_rdata
output_dir <- args$output_dir
name <- args$name

annotation_ucsc <- args$ucsc
annotation_file <- args$annotation_file
ucsc_genome <- args$ucscgenome
species <- args$species

# 1 load Rdata
load(ocr_rdata_file)
load(tss_rdata_file)

# 2 create txdb object
writeLines("\nCreating txdb ...")
genome_info <- GenomeInfoDb::keepStandardChromosomes(rtracklayer::SeqinfoForUCSCGenome(ucsc_genome), # nolint: line_length_linter.
                                                     species = species)

if (annotation_ucsc == "TRUE") {
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  writeLines("\nusing TxDb.Hsapiens.UCSC.hg38.knownGene")
} else if (annotation_ucsc == "FALSE") {
  txdb <- importGFFAnnotation(annotation_file, genome_info, species = species)
  writeLines("\nusing gencode annotation")
} else {
  print("need to provide annotation option")
}


# 3 annotation
writeLines("\nAnnotation ...")

writeLines("\nocr_train_gr Annotation")
ocr_train_gr <- CAGEfightR::assignTxType(ocr_train_gr, txModels = txdb)

writeLines("\nocr_train_aug_gr Annotation")
ocr_train_aug_gr <- CAGEfightR::assignTxType(ocr_train_aug_gr, txModels = txdb)

writeLines("\nocr_test_gr Annotation")
ocr_test_gr <- CAGEfightR::assignTxType(ocr_test_gr, txModels = txdb)

writeLines("\nocr_test_aug_gr Annotation")
ocr_test_aug_gr <- CAGEfightR::assignTxType(ocr_test_aug_gr, txModels = txdb)

writeLines("\nocrlike_neg_train_gr Annotation")
ocrlike_neg_train_gr <- CAGEfightR::assignTxType(ocrlike_neg_train_gr, txModels = txdb) # nolint: line_length_linter.

writeLines("\nocrlike_neg_train_aug_gr Annotation")
ocrlike_neg_train_aug_gr <- CAGEfightR::assignTxType(ocrlike_neg_train_aug_gr, txModels = txdb) # nolint: line_length_linter.

writeLines("\nocrlike_neg_test_gr Annotation")
ocrlike_neg_test_gr <- CAGEfightR::assignTxType(ocrlike_neg_test_gr, txModels = txdb) # nolint: line_length_linter.

writeLines("\nocrlike_neg_test_aug_gr Annotation")
ocrlike_neg_test_aug_gr <- CAGEfightR::assignTxType(ocrlike_neg_test_aug_gr, txModels = txdb) # nolint: line_length_linter.

writeLines("\nnearpos_neg_train_aug_gr Annotation")
nearpos_neg_train_aug_gr <- CAGEfightR::assignTxType(nearpos_neg_train_aug_gr, txModels = txdb) # nolint: line_length_linter.

writeLines("\nnearpos_neg_test_aug_gr Annotation")
nearpos_neg_test_aug_gr <- CAGEfightR::assignTxType(nearpos_neg_test_aug_gr, txModels = txdb) # nolint: line_length_linter.


# 4 save
#writeLines("\nSaving..")
#save(list = c("ocr_train_gr",
#              "ocr_train_aug_gr",
#              "ocr_test_gr",
#              "ocrlike_neg_train_gr",
#              "ocrlike_neg_train_aug_gr",
#              "ocrlike_neg_test_gr"),
#     file = paste0("4_", name, "_ocr_ocrlikeneg_augmentation_with_annotation.RData")) # nolint: line_length_linter.

# 4 save
writeLines("\nSaving..")
save(list = c("ocr_train_gr", "ocr_train_aug_gr",
              "ocr_test_gr", "ocr_test_aug_gr",
              "ocrlike_neg_train_gr", "ocrlike_neg_train_aug_gr",
              "ocrlike_neg_test_gr", "ocrlike_neg_test_aug_gr",
              "nearpos_neg_train_aug_gr", "nearpos_neg_test_aug_gr"),
     file = paste0("4_", name, "_augmentation_with_annotation.RData"))