# Load necessary libraries
library(CAGEfightR)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicRanges)
library(rtracklayer)
library(dplyr)

# Set input and output directories
input_dir <- "/Users/natsudanav/Desktop/PRIMEloci/evaluation/model_evaluation/data/GM12878_wt10M_profiles_te/metadata/"
output_dir <- "/Users/natsudanav/Desktop/PRIMEloci/evaluation/model_evaluation/data/GM12878_wt10M_profiles_te/metadata/out"

# Ensure output directory exists
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Load BSgenome for hg38
hg38 <- BSgenome.Hsapiens.UCSC.hg38

# Load UCSC TxDb annotations
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
#tx_gr <- genes(txdb)  # Extract gene annotations as GRanges

# Get the list of all CSV files in the input directory
files <- list.files(input_dir, pattern = "\\.csv$", full.names = TRUE)

for (file in files) {
  # Step 1: Read the CSV file
  df <- read.csv(file)

  # Step 2: Convert to GRanges object
  gr <- GRanges(
    seqnames = df$seqnames,
    ranges = IRanges(start = df$start, end = df$end),
    strand = df$strand
  )

  # Step 3: Annotate using CAGEfightR
  annotated_gr <- CAGEfightR::assignTxType(gr, txdb)

  # Step 4: Add annotation back to the data frame
  df$annotation <- annotated_gr$txType

  # Step 5: Save the annotated data frame
  output_file <- file.path(output_dir, paste0(basename(file)))
  write.csv(df, output_file, row.names = FALSE)

  cat("Annotated file saved to:", output_file, "\n")
}

cat("All files processed and saved in:", output_dir, "\n")
