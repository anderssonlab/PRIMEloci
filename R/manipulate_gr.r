#' Convert Strand Information to No Strand for GRanges Object
#'
#' This function takes a `GRanges` object
#' and sets the strand information to `"*"`,
#' indicating no strand specificity.
#'
#' @param region_gr A `GRanges` object containing genomic ranges.
#'
#' @return A `GRanges` object with the strand information set to `"*"`.
#'
#' @details The function modifies the strand information
#' of the input `GRanges` object, setting it to `"*"`.
#' This is useful in situations where strand information
#' is irrelevant or not needed.
#'
#' @importFrom GenomicRanges strand
#' @export
convert_strand_to_nostrand_gr <- function(region_gr) {
  GenomicRanges::strand(region_gr) <- "*"
  return(region_gr)
}


#' Remove Metadata and Duplicate Genomic Ranges
#'
#' This function takes a `GRanges` object, removes all metadata,
#' and then eliminates duplicate genomic ranges based on sequence names,
#' start and end positions, and strand information.
#'
#' @param gr A `GRanges` object containing genomic ranges
#' with or without metadata.
#'
#' @return A `GRanges` object without any metadata and
#' without duplicate genomic ranges.
#'
#' @details The function will strip the input `GRanges` object
#' of all metadata columns and then identify
#' and remove duplicate genomic ranges.
#' Only the `seqnames`, `ranges`, and `strand` information
#' will be considered for identifying duplicates.
#'
#' @examples
#' # Create a GRanges object with metadata
#' gr <- GenomicRanges::GRanges(
#'     seqnames = c("chr1", "chr1", "chr2", "chr2", "chr3", "chr3"),
#'     ranges = IRanges::IRanges(start = c(100, 100, 200, 250, 300, 300),
#'                               end = c(150, 150, 250, 250, 350, 350)),
#'     strand = c("+", "+", "-", "-", "+", "+"),
#'     score = c(5.0, 5.0, 4.0, 4.0, 3.0, 3.0)  # Example metadata
#' )
#'
#' # Remove metadata and duplicate genomic ranges
#' unique_gr <- remove_metadata_and_duplicates(gr)
#'
#' # Inspect the result
#' unique_gr
#'
#' @importFrom GenomicRanges GRanges seqnames ranges strand start end duplicated
#' @importFrom IRanges IRanges
#' @export
remove_metadata_and_duplicates <- function(gr) {
  # Remove metadata columns by creating a new GRanges object without metadata
  gr_no_metadata <- GenomicRanges::GRanges(
    seqnames = GenomicRanges::seqnames(gr),
    ranges = IRanges::IRanges(start = GenomicRanges::start(gr), 
                              end = GenomicRanges::end(gr)),
    strand = GenomicRanges::strand(gr)
  )

  # Identify duplicated ranges based on seqnames, ranges, and strand
  duplicated_indices <- GenomicRanges::duplicated(gr_no_metadata)

  # Subset the GRanges object to keep only unique ranges
  unique_gr <- gr_no_metadata[!duplicated_indices]

  return(unique_gr)
}
