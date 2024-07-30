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


## 0. Setting

# Retrieve environment variables
ext_dis <- as.numeric(Sys.getenv("EXT_DIS"))

dir_resources <- Sys.getenv("DIR_RESOURCES")
dir_results <- Sys.getenv("DIR_RESULTS")

infile_epc <- Sys.getenv("INFILE_EPC")

outfile_pos_gr <- Sys.getenv("OUTFILE_POS_GR")
outfile_pos_csv <- Sys.getenv("OUTFILE_POS_CSV")

outfile_neg_notsig_gr <- Sys.getenv("OUTFILE_NEG_NOTSIG_GR")
outfile_neg_notsig_csv <- Sys.getenv("OUTFILE_NEG_NOTSIG_CSV")

# outfile_plot_txtype_percentage <- Sys.getenv("OUTFILE_PLOT_TXTYPE_PERCENTAGE") # nolint

# Genome information and standard chromosomes from UCSC
species <- Sys.getenv("SPECIES")
genome_info <- GenomeInfoDb::keepStandardChromosomes(
  rtracklayer::SeqinfoForUCSCGenome(Sys.getenv("GENOME_INFO")),
  species = species
)
std_chr <- GenomicRanges::seqnames(genome_info)

# Dynamically load TxDb package and create the TxDb object
txdb_package <- Sys.getenv("TXDB_PACKAGE")

if (!requireNamespace(txdb_package, quietly = TRUE)) {
  stop(paste("Package", txdb_package, "not installed"))
}

do.call("library", list(txdb_package))

txdb <- get(txdb_package, envir = asNamespace(txdb_package))



## 1. Process

# Import file, mutate (start, end), and keep std chr
epc_raw <- read_tsv(file.path(dir_resources, infile_epc),
                    show_col_types = FALSE)

# Select specific columns, and filter out NA for "Regulated" and "Significant"
epc_df <- epc_raw %>%
  data.frame() %>%
  dplyr::select(
    c("chrom", "chromStart", "chromEnd", "measuredGeneSymbol",
      "Regulated", "EffectSize", "Significant", "pValueAdjusted")
  ) %>%
  dplyr::filter(.data$Regulated != "NA", .data$Significant != "NA")

# Change column names to match functions
colnames(epc_df) <- c("seqnames", "start", "end", "measuredGeneSymbol",
                      "Regulated", "EffectSize", "Significant",
                      "pValueAdjusted")

# Keep only standard chromosomes
epc_df <- keep_given_chromosomes_df(epc_df, std_chr)

# Categorize data based on regulated and significant
pos_df <- epc_df %>% dplyr::filter(.data$Regulated == TRUE)
neg_notsig_df <- epc_df %>%
  dplyr::filter(.data$Regulated == FALSE, .data$Significant == FALSE)
neg_sigfcpos_df <- epc_df %>%
  dplyr::filter(.data$Regulated == FALSE,
                .data$Significant == TRUE,
                .data$EffectSize >= 0)

# Extend from the midpoint between start and end to the setting distance
# Grouping the 'exactly' same identified enhancers and convert to GRanges
pos_gr <- pos_df %>%
  extend_from_center_thick_df(., ext_dis) %>%
  group_exact_positions_df() %>%
  GRanges()
neg_notsig_gr <- neg_notsig_df %>%
  extend_from_center_thick_df(., ext_dis) %>%
  group_exact_positions_df() %>%
  GRanges()
neg_sigfcpos_gr <- neg_sigfcpos_df %>%
  extend_from_center_thick_df(., ext_dis) %>%
  group_exact_positions_df() %>%
  GRanges()

# Select non-overlap identified enhancers in negative groups
# that are also in positive group and negative significant group
# prioritization: pos > neg_sigfcpos > neg_notsig
# negative significant group is not considered in the following steps
neg_notsig_gr <- IRanges::subsetByOverlaps(neg_notsig_gr,
                                           c(pos_gr, neg_sigfcpos_gr),
                                           type = "any",
                                           invert = TRUE)

# Annotation
# Default was set (tssUpstream = 100, tssDownstream = 100)
pos_gr <- CAGEfightR::assignTxType(pos_gr, txModels = txdb)
neg_notsig_gr <- CAGEfightR::assignTxType(neg_notsig_gr, txModels = txdb)

# Save the GRanges object
saveRDS(pos_gr, file.path(dir_results, outfile_pos_gr))
write.csv(data.frame(pos_gr), file.path(dir_results, outfile_pos_csv))

saveRDS(neg_notsig_gr, file.path(dir_results, outfile_neg_notsig_gr))
write.csv(data.frame(neg_notsig_gr),
          file.path(dir_results, outfile_neg_notsig_csv))