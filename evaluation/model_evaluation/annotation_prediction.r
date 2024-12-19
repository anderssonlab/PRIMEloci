# Load necessary libraries
library(CAGEfightR)
library(BSgenome.Hsapiens.UCSC.hg38)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(GenomicRanges)
library(rtracklayer)

# File paths (adjust paths accordingly)
input_dir <- "/Users/natsudanav/Desktop/PRIMEloci/evaluation/model_evaluation/data/GM12878_wt10M_profiles_te/predictions/"
output_dir <- "/Users/natsudanav/Desktop/PRIMEloci/evaluation/model_evaluation/data/GM12878_wt10M_profiles_te/predictions/out"


# Ensure output directory exists
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Load UCSC hg38 annotations
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# Get the list of all BED files in the input directory
bed_files <- list.files(input_dir, pattern = "\\.bed$", full.names = TRUE)

# Loop through each BED file in the directory
for (file in bed_files) {
  # Step 1: Load the BED file as a data.frame
  bed_df <- read.table(file, header = TRUE, stringsAsFactors = FALSE)
  
  # Step 2: Convert to GRanges object with adjusted start (BED is 0-based)
  gr <- GRanges(
    seqnames = bed_df$chrom,
    ranges = IRanges(start = bed_df$chromStart + 1,  # Convert to 1-based
                     end = bed_df$chromEnd),
    strand = bed_df$strand
  )
  
  # Step 3: Annotate using CAGEfightR
  annotated_gr <- CAGEfightR::assignTxType(gr, txdb)
  
  # Step 4: Add annotations back to the original BED file
  annotated_df <- data.frame(
    chrom = as.character(seqnames(annotated_gr)),
    chromStart = start(annotated_gr) - 1,  # Convert back to 0-based
    chromEnd = end(annotated_gr),
    name = bed_df$name,
    score = bed_df$score,
    strand = as.character(strand(annotated_gr)),
    sum_count = bed_df$sum_count,
    annotation = annotated_gr$txType  # Add annotation
  )
  
  # Step 5: Save the annotated data frame with column names
  output_file <- file.path(output_dir, paste0(basename(file)))
  write.table(annotated_df, file = output_file, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  
  cat("Annotated BED file saved to:", output_file, "\n")
}

cat("All BED files processed and saved in:", output_dir, "\n")

