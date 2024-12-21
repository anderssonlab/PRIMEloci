library(GenomicRanges)

# Define the file path
bed_file <- "/Users/natsudanav/Documents/PRIMEloci_pred_0_75_FANTOM5_rmSingletons_combined_coreovlwith-d.bed"
name_outfile <- "FANTOM5_rmSingletons_0.75_gr"
# Read the BED file
bed_data <- read.delim(bed_file, header = TRUE, stringsAsFactors = FALSE)

# BED files are 0-based and half-open, so adjust to 1-based, fully closed
bed_data$thickStart <- bed_data$thickStart + 1  # Convert to 1-based indexing

# Extend thickStart to create a fixed 401 bp region
adjusted_start <- pmax(0, bed_data$thickStart - 200) # Start is 200 bp upstream
adjusted_end <- bed_data$thickStart + 200            # End is 200 bp downstream

# Create GRanges object
gr <- GRanges(
  seqnames = bed_data$chrom,
  ranges = IRanges(
    start = adjusted_start,
    end = adjusted_end,
    names = bed_data$name
  ),
  strand = bed_data$strand,
  thick = IRanges(start = bed_data$thickStart, end = bed_data$thickStart),
  max_score = bed_data$max_score
)

## Combine into a GRangesList with a named element
#grl <- GRangesList(FANTOM5_rmSingletons_0.75 = gr)

## Save the GRangesList object as .rds
#saveRDS(grl, file = paste0(name_outfile, ".rds"))
saveRDS(gr, file = paste0(name_outfile, ".rds"))
