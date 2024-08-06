#' Generate Random Integers with Normal Distribution
#'
#' This function generates a specified number of random integers from a normal
#' distribution, truncated to lie within the range -30 to 30.
#'
#' @param number An integer specifying the number of random integers to generate
#' @return A vector of random integers.
#' @export
#' @examples
#' rnorm_int(5)
rnorm_int <- function(number) {
  ch <- 0
  ls <- c()

  while (ch < number) {
    cntrl_num <- number - ch
    nn <- round(stats::rnorm(cntrl_num, mean = 0, sd = 10), 0)
    nn <- nn[nn >= -30 & nn <= 30]
    ls <- c(ls, nn)
    ch <- length(ls)
  }

  return(ls[1:number])
}

#' Generate Random Numbers with Uniform and Normal Distributions
#'
#' This function generates random integers from both uniform and normal
#' distributions and plots their histograms.
#'
#' @param n An integer specifying the number of random numbers to generate.
#' @return A list containing two vectors: "unif" for uniform distribution and
#' "norm" for normal distribution.
#' @export
#' @examples
#' get_random_number(10)
get_random_number <- function(n) {
  unif_int <- round(runif(n = n, min = -30, max = 30), 0)

  norm_int <- rnorm_int(n)

  return(list("unif" = unif_int, "norm" = norm_int))
}

#' Data Augmentation with Uniform Distribution
#'
#' This function applies uniform distribution shifts to genomic ranges in the
#' given GRanges object and filters the results to avoid overlaps.
#'
#' @param ocr A GRanges object with a 'thick' column of IRanges.
#' @return A GRanges object with uniformly shifted ranges.
#' @export
#' @importFrom IRanges shift subsetByOverlaps
#' @examples
#' # ocr is a GRanges object with required data
#' augmented_ocr <- augmentation_unif(ocr)
augmentation_unif <- function(ocr) {
  # Assert that the 'thick' column is present and is an IRanges object
  assertthat::assert_that("thick" %in% names(S4Vectors::mcols(ocr)),
                          msg = "The 'thick' column is not present in the GRanges object.") # nolint: line_length_linter.
  assertthat::assert_that(is(ocr$thick, "IRanges"),
                          msg = "The 'thick' column is not an IRanges object.")

  num <- length(ocr)
  set <- get_random_number(num)

  # ocr_unif
  ocr_unif <- ocr
  start(ocr_unif) <- start(ocr_unif) + set$unif
  end(ocr_unif) <- end(ocr_unif) + set$unif
  ocr_unif$thick <- IRanges::shift(ocr_unif$thick, shift = set$unif)
  ocr_unif <- IRanges::subsetByOverlaps(ocr_unif,
                                        ocr,
                                        type = "equal",
                                        invert = TRUE)

  return(ocr_unif)
}

#' Data Augmentation with Normal Distribution
#'
#' This function applies normal distribution shifts to genomic ranges in the
#' given GRanges object and filters the results to avoid overlaps.
#'
#' @param ocr A GRanges object with a 'thick' column of IRanges.
#' @return A GRanges object with normally shifted ranges.
#' @export
#' @importFrom IRanges shift subsetByOverlaps
#' @examples
#' # ocr is a GRanges object with required data
#' augmented_ocr <- augmentation_norm(ocr)
augmentation_norm <- function(ocr) {
  # Assert that the 'thick' column is present and is an IRanges object
  assertthat::assert_that("thick" %in% names(S4Vectors::mcols(ocr)),
                          msg = "The 'thick' column is not present in the GRanges object.") # nolint: line_length_linter.
  assertthat::assert_that(is(ocr$thick, "IRanges"),
                          msg = "The 'thick' column is not an IRanges object.")

  num <- length(ocr)
  set <- get_random_number(num)

  # ocr_norm
  ocr_norm <- ocr
  start(ocr_norm) <- start(ocr_norm) + set$norm
  end(ocr_norm) <- end(ocr_norm) + set$norm
  ocr_norm$thick <- IRanges::shift(ocr_norm$thick, shift = set$norm)
  ocr_norm <- IRanges::subsetByOverlaps(ocr_norm,
                                        ocr,
                                        type = "equal",
                                        invert = TRUE)

  return(ocr_norm)
}
