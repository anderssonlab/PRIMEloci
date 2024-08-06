#!/usr/bin/env Rscript

writeLines("\n### Running _filter_bed_to_reduce_gr.r ###")

writeLines("\n# Importing R libraries..")
suppressPackageStartupMessages({
  library(assertthat)
  library(data.table)
  library(argparse)
  library(foreach)
  library(doParallel)

  library(GenomicRanges)
  library(rtracklayer)
  library(PRIMEloci)
})

# Define the command line arguments
parser <- ArgumentParser()
parser$add_argument("-i", "--input_bed", type = "character",
                    help = "Input BED file")
parser$add_argument("-o", "--output_dir", type = "character",
                    default = FALSE,
                    help = "Output directory, default will be the same as the input file") # nolint: line_length_linter.
args <- parser$parse_args()

# Load the BED file into a data.table
bed_file <- fread(args$input_bed, header = TRUE)

required_cols <- c("chrom", "chromStart", "chromEnd", "strand", "score")
assert_that(all(required_cols %in% colnames(bed_file)),
            msg = "The BED file must contain 'chrom', 'chromStart', 'chromEnd', 'strand', and 'score' columns.") # nolint: line_length_linter.


# Create a GRanges object from the data.table, including all columns as metadata
# Create a GRanges object from the data.table
gr <- GRanges(seqnames = bed_file$chrom,
              ranges = IRanges(start = bed_file$chromStart,
                               end = bed_file$chromEnd),
              strand = bed_file$strand)

# Add the rest of the columns as metadata columns
mcols(gr) <- bed_file[, !c("chrom",
                           "chromStart",
                           "chromEnd",
                           "strand"),
                      with = FALSE]

writeLines("Starting to get highest non-overlapping ranges...")
writeLines("It may take a while depending on the size of the input file...")
start_time <- Sys.time()
print(start_time)

# Get the highest non-overlapping ranges
selected_gr <- get_highest_non_overlap(gr)

print(selected_gr)
end_time <- Sys.time()
print(end_time - start_time)

# Determine the output file name
if (args$output_dir != FALSE) {
  output_dir <- args$output_dir
} else {
  output_dir <- dirname(args$input_bed)
}

input_basename <- tools::file_path_sans_ext(basename(args$input_bed))
output_rds <- file.path(output_dir, paste0(input_basename, "_reduced.rds"))

# Save the selected GRanges object to an RDS file
saveRDS(selected_gr, file = output_rds)
cat("Reduced GRanges object saved to", output_rds, "\n")

# Export the selected GRanges to a BED file with required columns and metadata
bed_df <- as.data.frame(selected_gr)
bed_df <- bed_df[, c("seqnames",
                     "start",
                     "end",
                     "width",
                     "strand", colnames(mcols(selected_gr)))]
setnames(bed_df,
         c("seqnames", "start", "end", "strand"),
         c("chrom", "chromStart", "chromEnd", "strand"))
bed_df <- bed_df[, colnames(bed_file)]
data.table::setorder(bed_df, chrom, chromStart)

# Write to BED file
output_bed <- file.path(output_dir, paste0(input_basename, "_reduced.bed"))
data.table::fwrite(bed_df,
                   file = output_bed,
                   sep = "\t",
                   quote = FALSE,
                   col.names = TRUE)

cat("Reduced GRanges object saved to", output_bed, "\n")
writeLines("Done!")
