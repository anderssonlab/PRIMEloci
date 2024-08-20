#' Load a BED file and validate its columns
#'
#' This function reads a BED file into a `data.table` and checks that it contains
#' the required columns: 'chrom', 'chromStart', 'chromEnd', 'strand', and 'score'.
#'
#' @param input_bed Character. The path to the input BED file.
#'
#' @return A `data.table` containing the BED file data.
#'
#' @importFrom data.table fread
#' @importFrom assertthat assert_that
#' @export
load_bed_file <- function(input_bed) {
  bed_file <- read.table(input_bed, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  required_cols <- c("chrom", "chromStart", "chromEnd", "strand", "score")
  assertthat::assert_that(all(required_cols %in% colnames(bed_file)),
                          msg = "The BED file must contain 'chrom', 'chromStart', 'chromEnd', 'strand', and 'score' columns.")
  return(bed_file)
}

#' Create a GRanges object from a BED data.table
#'
#' This function converts a `data.table` containing BED file data into a `GRanges` object.
#' The required columns are extracted and used to define the `GRanges` object, and the remaining
#' columns are added as metadata.
#'
#' @param bed_file A `data.table` containing the BED file data.
#'
#' @return A `GRanges` object.
#'
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors mcols
#' @export
create_granges_from_bed <- function(bed_file) {
  gr <- GenomicRanges::GRanges(seqnames = bed_file$chrom,
                               ranges = IRanges::IRanges(start = bed_file$chromStart + 1,
                                                         end = bed_file$chromEnd),
                               strand = bed_file$strand)
  S4Vectors::mcols(gr) <- bed_file[, !(names(bed_file) %in% c("chrom", "chromStart", "chromEnd", "strand"))]
  return(gr)
}

#' Save a GRanges object to RDS and BED formats
#'
#' This function saves a `GRanges` object as an RDS file and converts it to a `data.frame`
#' for saving as a BED file. The BED file is saved with adjusted columns to match the original format.
#'
#' @param gr A `GRanges` object to be saved.
#' @param output_dir Character. The directory where the output files will be saved.
#' @param input_basename Character. The base name for the output files.
#'
#' @return None. The function saves the output files to the specified directory.
#'
#' @importFrom tools file_path_sans_ext
#' @importFrom data.table fwrite
#' @importFrom S4Vectors mcols
#' @export
save_granges_to_bed <- function(gr, output_dir, input_basename, bed_file) {
  # Save the selected GRanges object to an RDS file
  output_rds <- file.path(output_dir, paste0(input_basename, "_reduced.rds"))
  saveRDS(gr, file = output_rds)
  cat("Reduced GRanges object saved to", output_rds, "\n")

  # Convert GRanges object to a data frame
  bed_df <- as.data.frame(gr)

  # Select and rearrange the columns
  bed_df <- bed_df[, c("seqnames", "start", "end", "width", "strand", colnames(S4Vectors::mcols(gr)))]

  # Adjust the 'start' column to be 0-based by subtracting 1
  bed_df$start <- bed_df$start - 1

  # Rename the columns to match BED file format
  data.table::setnames(bed_df, c("seqnames", "start", "end", "strand"), c("chrom", "chromStart", "chromEnd", "strand"))

  # Rearrange the columns to match the original BED file format
  bed_df <- bed_df[, colnames(bed_file)]

  # Sort the data frame by chrom and chromStart
  bed_df <- bed_df[order(bed_df$chrom, bed_df$chromStart), ]

  # Write to BED file
  output_bed <- file.path(output_dir, paste0(input_basename, "_reduced.bed"))
  data.table::fwrite(bed_df, file = output_bed, sep = "\t", quote = FALSE, col.names = TRUE)

  cat("Reduced GRanges object saved to", output_bed, "\n")
}