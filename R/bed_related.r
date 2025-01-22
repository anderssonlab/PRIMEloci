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