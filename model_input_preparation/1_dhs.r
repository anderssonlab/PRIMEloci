#!/usr/bin/env Rscript

writeLines("\n #### 1 PREPARING DHS DATA ####")

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


### ARGPARSE
parser <- ArgumentParser()

parser$add_argument("-i", "--input_file",
                    default = "filename",
                    help = "dhs input files path")
parser$add_argument("-n", "--name", default = "name", help = "Name")
parser$add_argument("-o", "--output_dir",
                    default = "./",
                    help = "output directory")
parser$add_argument("-b", "--blacklist_file",
                    default = "../data/resources/hg38-blacklist.v2.bed",
                    help = "blacklist file path")
parser$add_argument("-d", "--distance",
                    default = 200, type = "integer",
                    help = "distance from the dhs peak (bps)")
parser$add_argument("-t", "--test_chr",
                    default = "chr2,chr3,chr4",
                    help = "chromosome(s) that you would like to keep for the testing set (no space between , )")  # nolint: line_length_linter.
parser$add_argument("-g", "--ucscgenome",
                    default = "hg38",
                    help = "UCSC genome")
parser$add_argument("-s", "--species",
                    default = "Homo sapiens",
                    help = "species")

args <- parser$parse_args()


# General
dist <- args$distance
test_chr <- unlist(strsplit(args$test_chr, ","))

# Input and output
input_file <- args$input_file
blacklist_file <- args$blacklist_file
output_dir <- args$output_dir
name <- args$nam
ucsc_genome <- args$ucscgenome
species <- args$species


# 1. Read in data
writeLines("\nReading in data and filtering blacklist..")

dhs_data_path <- list()
dhs_data_path[[name]] <- input_file

readin_dhs_data <- import_dhs_encode(dhs_data_path)

# Load blacklist and remove overlapping dhs regions
# https://github.com/Boyle-Lab/Blacklist
blacklist <- import(blacklist_file, format = "bed")

# Extend +- from midposition, and keep only standard chromosomes
genome_info <- GenomeInfoDb::keepStandardChromosomes(rtracklayer::SeqinfoForUCSCGenome(ucsc_genome), # nolint: line_length_linter.
                                                     species = species)
genome_info <- GenomeInfoDb::dropSeqlevels(genome_info, "chrM")
standard_chromosomes <- GenomeInfoDb::seqnames(genome_info)
print("Standard chromosomes:")
print(standard_chromosomes)

dhs_data <- lapply(readin_dhs_data, function(df) {
  df %>%
    extend_from_center_thick_df(.,
                                dist,
                                keep_same_length = TRUE,
                                group_exact_positions = TRUE) %>%
    dplyr::filter(seqnames %in% standard_chromosomes)
})

dhs_data <- purrr::map2_df(dhs_data,
                           names(dhs_data),
                           ~dplyr::mutate(.x, lib_name = .y))


# 2. Training: filter out the test chromosomes, and remove blacklist
ocr_train_gr <- dhs_data %>%
  dplyr::filter(!(seqnames %in% test_chr)) %>%
  GRanges() %>%
  GRanges(thick = IRanges(start = .$thick, end = .$thick))
ocr_train_gr <- filter_blacklist(ocr_train_gr, blacklist)
names(ocr_train_gr) <- paste0(seqnames(ocr_train_gr), ":",
                              start(ocr_train_gr), "-",
                              end(ocr_train_gr), ";",
                              strand(ocr_train_gr))

# 3. Testing: select the test chromosomes, and remove blacklist
ocr_test_gr <- dhs_data %>%
  dplyr::filter(seqnames %in% test_chr) %>%
  GRanges(thick = IRanges(start = .$thick, end = .$thick))
ocr_test_gr <- filter_blacklist(ocr_test_gr, blacklist)
names(ocr_test_gr) <- paste0(seqnames(ocr_test_gr), ":",
                             start(ocr_test_gr), "-",
                             end(ocr_test_gr), ";",
                             strand(ocr_test_gr))

writeLines("\nReport..")
print(paste0("Input dhs file: ", input_file))
print(paste0("Blacklist file: ", blacklist_file))
print(paste0("Output: ", "1.", name, ".ocr.RData"))
print(paste0("Distance from the dhs peak (bps): ", dist))
print(paste0("OCR train: ", length(ocr_train_gr)))
print(paste0("OCR test: ", length(ocr_test_gr)))

# 4. Save
writeLines("\nSaving..")
save(list = c("ocr_train_gr", "ocr_test_gr"),
     file = paste0(output_dir, "1_", name, "_ocr.RData"))
