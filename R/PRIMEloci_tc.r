#' Get tag clusters and extend from thick positions
#'
#' This function calculates tag clusters from
#' a RangedSummarizedExperiment object \code{ctss_rse}
#' and extends them by a specified distance \code{ext_dis}
#' around the thick positions.
#'
#' @param ctss_rse A RangedSummarizedExperiment object containing CAGE data.
#' @param ext_dis An integer specifying the distance
#' to extend around thick positions (default is 200).
#'
#' @return A GenomicRanges::GRangesList object where each element corresponds to
#'         a tag cluster extended by \code{ext_dis} around the thick positions.
#'
#' @details
#' The function iterates over columns of \code{ctss_rse},
#' calculates pooled counts,
#' subsets clusters based on a score threshold (> 0),
#' and clusters them unidirectionally.
#' It then extends each cluster by \code{ext_dis} base pairs
#' around the thick positions
#' and stores the results in a GRangesList object.
#'
#' @examples
#' # Example usage with a RangedSummarizedExperiment object
#' # ctss_rse <- ...  # Load or create your RangedSummarizedExperiment object
#' # ext_dis <- 200   # Define your extension distance
#' # result <- get_tagclusters_and_extend_fromthick(ctss_rse, ext_dis)
#'
#' @import SummarizedExperiment
#' @import GenomicRanges
#' @import IRanges
#' @import CAGEfightR
#' @import assertthat
#'
#' @export
get_tcs_and_extend_fromthick <- function(ctss_rse, ext_dis = 200) {

  # Assert that ctss_rse is a RangedSummarizedExperiment object
  assertthat::assert_that(
    inherits(ctss_rse, "RangedSummarizedExperiment"),
    msg = "ctss_rse must be a RangedSummarizedExperiment object."
  )

  # Get column names
  col_ctss_rse <- colnames(ctss_rse)

  # Initialize GRangesList
  tc_grl <- GenomicRanges::GRangesList()

  # Loop over column names
  for (i in col_ctss_rse) {

    writeLines(paste("Processing: ", i, "\n"))

    # Extract data for current column
    ctss <- ctss_rse[, i]

    # Calculate pooled counts
    ctss <- CAGEfightR::calcPooled(ctss, inputAssay = "counts")

    # Subset by score > 0 (score column is created by calcPooled())
    ctss <- base::subset(ctss, score > 0)

    # Cluster unidirectionally
    object <- CAGEfightR::clusterUnidirectionally(ctss)

    # Create new ranges around thick positions
    new_ranges <- IRanges::IRanges(start = start(object$thick) - ext_dis,
                                   end = end(object$thick) + ext_dis)

    # Assign new ranges to object
    new_object <- object
    IRanges::ranges(new_object) <- new_ranges

    # Trim new object
    writeLines("Trimming out-of-bound ranges...")
    new_object <- GenomicRanges::trim(new_object)

    writeLines("Keep only prefered width...\n")
    len_vec <- ext_dis * 2 + 1

    new_object_widths <- GenomicRanges::width(new_object)
    new_object <- new_object[new_object_widths == len_vec]

    # Store in GRangesList
    tc_grl[[i]] <- new_object
  }

  return(tc_grl)
}


#' Validate a TC object
#'
#' This function ensures that a given `GenomicRanges::GRanges` or
#' `GenomicRanges::GRangesList` object meets the expected criteria.
#'
#' @param tc_object A `GenomicRanges::GRanges` or
#' `GenomicRanges::GRangesList` object.
#' @param ctss_rse The original CTSS dataset (for length and naming validation).
#' @param ext_dis An integer value for extension distance. Default is 200.
#'
#' @return TRUE if validation passes; otherwise, throws an error.
#'
#' #' @examples
#' # Validate an existing TC object:
#' validate_tc_object(tc_object, ctss_rse)
#'
#' @import GenomicRanges
#' @import assertthat
#'
#' @export
#'
validate_tc_object <- function(tc_object, ctss_rse, ext_dis = 200) {

  len_vec <- ext_dis * 2 + 1

  # Check if tc_object is a valid GRanges or GRangesList
  assertthat::assert_that(
    inherits(tc_object, "GenomicRanges::GRanges") || inherits(tc_object, "GenomicRanges::GRangesList"), # nolint: line_length_linter.
    msg = "tc_object must be a GenomicRanges::GRanges or GenomicRanges::GRangesList object" # nolint: line_length_linter.
  )

  if (inherits(tc_object, "GenomicRanges::GRanges")) {
    # Ensure all regions have the correct width
    assertthat::assert_that(
      all(GenomicRanges::width(tc_object) == len_vec),
      msg = paste("All regions in tc_object (GenomicRanges::GRanges) must have width", len_vec) # nolint: line_length_linter.
    )
  }

  if (inherits(tc_object, "GenomicRanges::GRangesList")) {
    # Ensure all GRanges in GRangesList have correct widths
    assertthat::assert_that(
      all(sapply(tc_object,
                 function(gr) all(GenomicRanges::width(gr) == len_vec))),
      msg = paste("All regions in each GRanges of tc_object (GenomicRanges::GRangesList) must have width", len_vec) # nolint: line_length_linter.
    )

    # Ensure the length of tc_object matches ctss_rse
    assertthat::assert_that(
      length(tc_object) == length(ctss_rse),
      msg = "tc_object (as GenomicRanges::GRangesList) must have the same length as ctss_rse" # nolint: line_length_linter.
    )

    # Ensure names match
    assertthat::assert_that(
      identical(names(tc_object), names(ctss_rse)),
      msg = "tc_object (as GenomicRanges::GRangesList) must have the same names as ctss_rse" # nolint: line_length_linter.
    )
  }

  message("TC object validation passed successfully.")
  return(TRUE)
}
