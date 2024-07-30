#!/usr/bin/env Rscript


# Load required libraries
library(GenomicRanges)
library(argparse)
library(data.table)

# Define the command line arguments
parser <- ArgumentParser()

parser$add_argument("-i", "--input_bed", type = "character",
                    help = "Input BED file")
parser$add_argument("-o", "--output_dir", type = "character",
                    default = "./", help = "Output directory")

args <- parser$parse_args()


# Load the BED file into a data.table
bed_file <- fread(args$input_bed, header = TRUE)

# Create a GRanges object from the data.table
gr <- GRanges(seqnames = bed_file$chrom,
              ranges = IRanges(start = bed_file$chromStart,
                               end = bed_file$chromEnd),
              strand = bed_file$strand,
              score = bed_file$score,
              name = bed_file$name)


# Find overlapping regions and reduce them,
# selecting the one with the highest score
reduced_gr <- reduce(gr, with.revmap = TRUE)
selected_gr <- unlist(GRangesList(lapply(reduced_gr$revmap, function(idx) {
  gr_subset <- gr[idx]
  gr_subset[which.max(gr_subset$score)]
})))


# Determine the output file name
input_basename <- tools::file_path_sans_ext(basename(args$input_bed))
output_rds <- file.path(args$output_dir, paste0(input_basename, "_reduced.rds"))

# Save the output to an RDS file
saveRDS(selected_gr, file = output_rds)

cat("Reduced GRanges object saved to", output_rds, "\n")
