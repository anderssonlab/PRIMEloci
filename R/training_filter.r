# source("../R/general.r")



#' Select OCRs Based on Specific Conditions for Each Library
#'
#' This function selects OCRs (open chromatin regions) from a CAGE dataset
#' based on specified conditions for each library. It filters the regions
#' based on the number of overlaps and scores, ensuring that selected OCRs
#' meet a minimum score threshold for one-base overlaps.
#'
#' @param cage_rse A `SummarizedExperiment` object containing the CAGE dataset.
#' @param ocr_gr A `GRanges` object representing open chromatin regions (OCRs).
#' @param coln An integer indicating the column index of the assay to use.
#' @param min_score_for_onebase_ocr A numeric value specifying the minimum
#'   score required for one-base overlaps. Defaults to 4.
#' @return A `GRanges` object containing the selected OCRs based on the
#'   specified conditions.
#' @importFrom GenomicRanges countOverlaps findOverlaps GRanges
#' @importFrom dplyr filter select arrange
#' @examples
#' # Example usage:
#' cage_rse <- ... # Define the SummarizedExperiment object
#' ocr_gr <- ... # Define the GRanges object for OCRs
#' selected_ocr <- select_ocr_on_condition_each_lib(cage_rse, ocr_gr, coln = 1)
#' @details
#' This function calls `cast_rse_to_granges` from the `general.r` script.
#' @export
select_ocr_on_condition_each_lib <- function(cage_rse,
                                             ocr_gr,
                                             coln,
                                             min_score_for_onebase_ocr = 4) {
  cage_gr <- cast_rse_to_granges(cage_rse, assay = "counts", coln_assay = coln)
  cage_filtered_gr <- cage_gr[cage_gr$score > 0]
  expr <- GenomicRanges::countOverlaps(ocr_gr, cage_filtered_gr)
  ocr_df <- data.frame(ocr_gr)
  ocr_df$num_overlaps <- expr

  sure_ocr_df <- dplyr::filter(ocr_df, num_overlaps > 1)

  onebase_ocr_df <- dplyr::filter(ocr_df, num_overlaps == 1)
  onebase_ocr_gr <- GenomicRanges::GRanges(onebase_ocr_df)
  onebase_ovl <- GenomicRanges::findOverlaps(onebase_ocr_gr, cage_filtered_gr)

  cage_hit_gr <- cage_filtered_gr[S4Vectors::subjectHits(onebase_ovl), ]
  onebase_ocr_df$onebase_score <- cage_hit_gr$score

  slt_onebase_ocr_df <- dplyr::filter(onebase_ocr_df, onebase_score >= min_score_for_onebase_ocr)
  slt_onebase_ocr_df <- dplyr::select(slt_onebase_ocr_df, -onebase_score)

  slt_ocr_df <- rbind(sure_ocr_df, slt_onebase_ocr_df)
  slt_ocr_df <- dplyr::arrange(slt_ocr_df, seqnames, start, end)

  return(GenomicRanges::GRanges(slt_ocr_df))
}


#' Select OCRs Based on Specific Conditions
#'
#' This function selects OCRs (open chromatin regions) from a CAGE dataset
#' across all libraries based on specified conditions. It uses the
#' `select_ocr_on_condition_each_lib` function for each library to ensure
#' consistency and accuracy in the selection process.
#'
#' @param cage_rse A `SummarizedExperiment` object containing the CAGE dataset.
#' @param ocr_gr A `GRanges` object representing open chromatin regions (OCRs).
#' @return A `GRangesList` object containing the selected OCRs for each library
#'   based on the specified conditions.
#' @importFrom GenomicRanges GRangesList
#' @examples
#' # Example usage:
#' cage_rse <- ... # Define the SummarizedExperiment object
#' ocr_gr <- ... # Define the GRanges object for OCRs
#' selected_ocr_grl <- select_ocr_on_condition(cage_rse, ocr_gr)
#' @details
#' This function iterates over all libraries in the CAGE dataset and applies
#' the selection criteria using `select_ocr_on_condition_each_lib`.
#' @export
select_ocr_on_condition <- function(cage_rse, ocr_gr) {
  slt_ocr_grl <- GenomicRanges::GRangesList()
  for (i in 1:ncol(cage_rse)) {
    slt_ocr_grl[[i]] <- select_ocr_on_condition_each_lib(cage_rse, ocr_gr, i)
  }
  return(slt_ocr_grl)
}
