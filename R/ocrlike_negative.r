# source("../R/profiles.r")
# source("../R/manipulate_gr.r")

# Character of negative
# 1. Singleton was removed before calling TC
# 2. TC won't overlap with the positive ocr
# 3. Extend DIST from THICK of TC
# It means that TC won't overlap with ocr, but extended TC might
# = we allow the neg to overlap with pos
# Note that TCs come from all selected CAGE.

#' Generate OCR-like Negative Samples
#'
#' This function generates OCR-like negative samples
#' by identifying transcription clusters (TCs)
#' from CAGE data that do not overlap with positive OCR regions.
#' The identified TCs are extended from their center points based on the "thick"
#' positions. The negative samples can overlap with positive regions after
#' extension, but the original TCs will not overlap with positive OCR regions.
#'
#' @param orig_ctss A CAGE dataset containing the original CTSS information.
#' @param ocr_gr A `GRanges` object representing positive OCR regions.
#' @param dis A numeric value specifying the distance to extend from the center
#'   of the TCs.
#' @return A `GRanges` object containing the OCR-like negative samples.
#' @importFrom CAGEfightR subsetBySupport calcPooled clusterUnidirectionally
#' @importFrom IRanges subsetByOverlaps
#' @importFrom GenomicRanges duplicated
#' @importFrom S4Vectors Rle
#' @examples
#' # Example usage:
#' # Assuming `orig_ctss` is a CAGE dataset and
#' `ocr_gr` is a GRanges object for OCR
#' negative_samples <- ocrlike_for_negative(orig_ctss, ocr_gr, dis = 100)
#' @export
ocrlike_for_negative <- function(orig_ctss, ocr_gr, dis) {
  supp_ctss <- CAGEfightR::subsetBySupport(orig_ctss,
                                           inputAssay = "counts",
                                           outputColumn = "support",
                                           unexpressed = 1,
                                           minSamples = 0)
  supp_ctss <- CAGEfightR::calcPooled(supp_ctss,
                                      inputAssay = "counts")
  tcs <- CAGEfightR::clusterUnidirectionally(supp_ctss)
  tcs_neg <- IRanges::subsetByOverlaps(tcs,
                                       ocr_gr,
                                       invert = TRUE)

  # Extend the negative (ocr-like) positions from the center of the thick position # nolint: line_length_linter.
  tcs_neg <- extend_from_center_thick_gr(tcs_neg, dis, keep_same_length = TRUE)
  mcols(tcs_neg) <- mcols(tcs_neg)[, "thick", drop = FALSE]

  # update row names to remove strand information
  strand(tcs_neg) <- S4Vectors::Rle(rep("*", length(tcs_neg)))

  # Remove duplicated genomic ranges
  duplicated_indices <- GenomicRanges::duplicated(tcs_neg)

  # Subset the GRanges object to keep only unique ranges
  print("Removing duplicated genomic ranges")
  unique_tcs_neg <- tcs_neg[!duplicated_indices]
  names(unique_tcs_neg) <- paste0(seqnames(unique_tcs_neg), ":",
                                  start(unique_tcs_neg), "-",
                                  end(unique_tcs_neg), ";",
                                  strand(unique_tcs_neg))

  # Report
  len_neg <- length(unique_tcs_neg)
  print(paste0("Number of ocr-like positions: ", len_neg))
  print(paste0("Ocr-like positions that (any) overlaps with the ocr: ", length(IRanges::subsetByOverlaps(unique_tcs_neg, ocr_gr, type = c("any")))))  # nolint: line_length_linter.
  print(paste0("Sanity check to make sure that no exact overlap of negative (ocr-like) to positive (ocr) (expect to be 0): ", length(IRanges::subsetByOverlaps(unique_tcs_neg, ocr_gr, type = c("equal"))))) # nolint: line_length_linter.

  return(unique_tcs_neg)

}