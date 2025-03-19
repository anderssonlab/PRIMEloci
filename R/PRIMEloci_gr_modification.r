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
#'
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
#'
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



#' Extend Genomic Ranges from the Center Based on Thickness
#'
#' This function extends genomic ranges from their center points, defined by
#' the "thick" positions in the metadata columns of a `GRanges` object. The
#' "thick" column must exist and be of class `IRanges`. The extension can
#' either maintain the original interval lengths or allow them to vary.
#'
#' @param gr A `GRanges` object containing genomic ranges with a "thick" column
#'   in the metadata.
#' @param dis A numeric value specifying the distance to extend from the center.
#' @param keep_same_length A logical value indicating whether to keep the
#'   intervals at the same length after extension. Defaults to TRUE.
#' @return A `GRanges` object with updated start and end positions, and an
#'   updated "thick" column if `keep_same_length` is TRUE.
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors mcols
#' @examples
#' # Assuming `gr` is a GRanges object with a 'thick' column of class IRanges
#' extended_gr <- extend_from_center_thick_gr(gr, dis = 100, keep_same_length = TRUE) # nolint: line_length_linter.
#'
extend_from_center_thick_gr <- function(gr,
                                        dis,
                                        keep_same_length = TRUE) {

  # Check if the thick column exists in the metadata
  if (!("thick" %in% names(S4Vectors::mcols(gr)))) {
    stop("The 'thick' column must exist in the metadata columns.")
  }

  # Check if the thick column is of class IRanges
  if (!inherits(S4Vectors::mcols(gr)$thick, "IRanges")) {
    stop("The 'thick' column must be of class IRanges.")
  }

  # If both checks pass, print a success message
  print("The 'thick' column exists and is of class IRanges.")

  start(gr) <- pmax(start(gr$thick) - dis, 1)

  if (keep_same_length) {
    end(gr) <- start(gr) + (dis * 2)
    mcols(gr)$thick <- IRanges::IRanges(start = start(gr) + dis,
                                        end = start(gr) + dis)
    print("If the extended start is out-of-range, the thick will be rewritten to the center of the extended negative range.") # nolint: line_length_linter.
  } else {
    end(gr) <- start(gr$thick) + dis
    print("The length of the intervals could be different. If you want to keep the same length, set 'keep_same_length = TRUE'.") # nolint: line_length_linter.
  }

  return(gr)
}



#' Remove 'chrM' from GRangesList
#'
#' This function removes all genomic ranges located on the chromosome 'chrM'
#' from each `GRanges` object within a `GRangesList`. It is useful for excluding
#' mitochondrial DNA, which is often represented by 'chrM', from analyses.
#'
#' @param grl A `GRangesList` object containing genomic ranges.
#' @return A `GRangesList` object with all ranges on 'chrM' removed.
#' @importFrom GenomicRanges GRangesList seqnames
#' @importFrom S4Vectors Rle
#' @examples
#' library(GenomicRanges)
#' gr1 <- GRanges(seqnames = c("chr1", "chrM", "chr2"),
#'                ranges = IRanges(start = 1:3, width = 3), strand = "+")
#' gr2 <- GRanges(seqnames = c("chrM", "chr2", "chr3"),
#'                ranges = IRanges(start = 4:6, width = 3), strand = "-")
#' grl <- GRangesList(gr1 = gr1, gr2 = gr2)
#'
#' # Remove 'chrM' ranges
#' modified_grl <- drop_chrM_from_grl(grl)
#'
drop_chrM_from_grl <- function(grl) { # nolint: object_name_linter.
  modified_grl <- lapply(grl, function(gr) {
    gr <- gr[seqnames(gr) != "chrM"]
    return(gr)
  })
  return(GenomicRanges::GRangesList(modified_grl))
}
