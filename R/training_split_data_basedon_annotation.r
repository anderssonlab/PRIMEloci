#' Subset Positive Annotations
#'
#' This function subsets positive annotations from a `GRangesList` object by 
#' removing entries with specific transcript types (`txType`), such as 'exon', 
#' 'CDS', and 'threeUTR'.
#'
#' @param pos_sltocr_grl A `GRangesList` object containing positive selected 
#'   open chromatin regions (OCRs) with a metadata column `txType` indicating 
#'   transcript types.
#' @return A `GRangesList` object with entries filtered to exclude those with 
#'   `txType` values 'exon', 'CDS', and 'threeUTR'.
#' @importFrom GenomicRanges GRangesList mcols
#' @importFrom S4Vectors subset
#' @importFrom assertthat assert_that
#' @examples
#' # Example usage:
#' # pos_sltocr_grl is a GRangesList object with a 'txType' metadata column
#' filtered_grl <- subset_annotation_positive(pos_sltocr_grl)
#' @export
subset_annotation_positive <- function(pos_sltocr_grl) {
  grl <- GenomicRanges::GRangesList()
  for (i in seq_along(pos_sltocr_grl)) {
    # Assert that txType is present in the GRanges object
    assertthat::assert_that("txType" %in% names(GenomicRanges::mcols(pos_sltocr_grl[[i]])),
                            msg = paste("txType is not present in the GRanges object at index", i))

    grl[[i]] <- S4Vectors::subset(pos_sltocr_grl[[i]], !(txType %in% c('exon', 'CDS', 'threeUTR')))
  }
  return(grl)
}



#' Select Negative Samples to Match Positive Samples
#'
#' This function selects a subset of negative samples to match the number of
#' positive samples in a given `GRangesList`. The selection is stratified by
#' chromosome (`seqnames`) and transcript type (`txType`) to ensure a balanced
#' representation.
#'
#' @param neg_sltocr_grl A `GRangesList` object containing negative selected
#'   open chromatin regions (OCRs).
#' @param match_pos_grl A `GRangesList` object containing positive selected
#'   OCRs to match in terms of sample size.
#' @return A `GRangesList` object containing a subset of the negative samples
#'   matching the size and distribution of the positive samples.
#' @importFrom GenomicRanges GRangesList GRanges
#' @importFrom dplyr group_by group_split sample_frac bind_rows
#' @importFrom purrr map
#' @examples
#' # Example usage:
#' # neg_sltocr_grl and match_pos_grl are GRangesList objects
#' selected_neg_grl <- select_negative_to_match_positive(neg_sltocr_grl, match_pos_grl)
#' @export
select_negative_to_match_positive <- function(neg_sltocr_grl, match_pos_grl) {
  grl <- GenomicRanges::GRangesList()

  for (i in seq_along(neg_sltocr_grl)) {
    new_df <- data.frame(neg_sltocr_grl[[i]])
    num_pos <- length(match_pos_grl[[i]])
    neg_ratio <- num_pos / nrow(new_df)

    split_data <- new_df %>%
      dplyr::group_by(seqnames, txType) %>%
      dplyr::group_split() %>%
      purrr::map(function(group_data) {
        sample_split <- dplyr::sample_frac(group_data, neg_ratio)
        list(negative = sample_split)
      })

    split_neg <- purrr::map(split_data, function(split_element) split_element$negative) %>% 
      dplyr::bind_rows()

    grl[[i]] <- GenomicRanges::GRanges(split_neg)
  }
  return(grl)
}


#' Split Data into Training and Validation Sets
#'
#' This function splits a given data frame into training and validation sets
#' based on a specified validation ratio. The split is stratified by
#' `seqnames` (chromosomes) and `txType` (transcript type).
#'
#' @param df A data frame containing the data to be split, with columns
#'   including `seqnames` and `txType`.
#' @param ratio_valid A numeric value specifying the proportion of the data
#'   to be used as the validation set. Defaults to 0.2.
#' @return A `GRangesList` object containing two `GRanges` objects: one for
#'   the training set and one for the validation set.
#' @importFrom dplyr group_by group_split sample_frac anti_join bind_rows arrange
#' @importFrom purrr map
#' @importFrom GenomicRanges GRanges GRangesList
#' @examples
#' df <- data.frame(seqnames = rep("chr1", 100), txType = rep("gene", 100), start = 1:100, end = 101:200)
#' split_sets <- split_train_valid(df, ratio_valid = 0.2)
#' @export
split_train_valid <- function(df, ratio_valid = 0.2) {
  # split based on chr and txType
  split_data <- df %>%
    dplyr::group_by(seqnames, txType) %>%
    dplyr::group_split() %>%
    purrr::map(function(group_data) {
      sample_split <- dplyr::sample_frac(group_data, ratio_valid)
      list(train = dplyr::anti_join(group_data, sample_split),
           valid = sample_split)
    })

  n <- length(split_data)
  print(n)
  valid_gr <- purrr::map(split_data[1:n], function(split_element) split_element$valid) %>%
    dplyr::bind_rows() %>%
    dplyr::arrange(seqnames, start, end) %>%
    GenomicRanges::GRanges()
  train_gr <- purrr::map(split_data[1:n], function(split_element) split_element$train) %>% 
    dplyr::bind_rows() %>%
    dplyr::arrange(seqnames, start, end) %>%
    GenomicRanges::GRanges()

  return(GenomicRanges::GRangesList(train_gr, valid_gr))
}


#' Wrap Up Data Splitting into Training and Validation Sets
#'
#' This function wraps up the process of splitting positive and negative
#' genomic ranges into training and validation sets. It ensures the number of
#' training and validation samples is balanced between positive and negative
#' sets.
#'
#' @param pos_grl A `GRangesList` object containing positive genomic ranges.
#' @param neg_grl A `GRangesList` object containing negative genomic ranges.
#' @return A list containing four `GRangesList` objects: `pos_tr_grl` for
#'   positive training sets, `pos_vl_grl` for positive validation sets,
#'   `neg_tr_grl` for negative training sets, and `neg_vl_grl` for negative
#'   validation sets.
#' @importFrom GenomicRanges GRangesList
#' @examples
#' pos_grl <- GenomicRanges::GRangesList()
#' neg_grl <- GenomicRanges::GRangesList()
#' result <- wrapup_split_train_valid(pos_grl, neg_grl)
#' @export
wrapup_split_train_valid <- function(pos_grl, neg_grl) {
  pos_tr_grl <- GenomicRanges::GRangesList()
  pos_vl_grl <- GenomicRanges::GRangesList()

  neg_tr_grl <- GenomicRanges::GRangesList()
  neg_vl_grl <- GenomicRanges::GRangesList()

  for (i in 1:length(pos_grl)) {
    pos_df <- data.frame(pos_grl[[i]])
    neg_df <- data.frame(neg_grl[[i]])

    pos_splt_grl <- split_train_valid(pos_df, ratio_valid = validation_ratio)
    neg_splt_grl <- split_train_valid(neg_df, ratio_valid = validation_ratio)

    final_num_tr <- min(length(pos_splt_grl[[1]]), length(neg_splt_grl[[1]]))
    final_num_vl <- min(length(pos_splt_grl[[2]]), length(neg_splt_grl[[2]]))

    pos_tr_grl[[i]] <- keep_equal_num(pos_splt_grl[[1]], final_num_tr)
    pos_vl_grl[[i]] <- keep_equal_num(pos_splt_grl[[2]], final_num_vl)

    neg_tr_grl[[i]] <- keep_equal_num(neg_splt_grl[[1]], final_num_tr)
    neg_vl_grl[[i]] <- keep_equal_num(neg_splt_grl[[2]], final_num_vl)
  }

  pos_tr_grl <- convert_to_thick_iranges(pos_tr_grl)
  pos_vl_grl <- convert_to_thick_iranges(pos_vl_grl)
  neg_tr_grl <- convert_to_thick_iranges(neg_tr_grl)
  neg_vl_grl <- convert_to_thick_iranges(neg_vl_grl)

  grl_list <- list(
    pos_tr_grl = pos_tr_grl,
    pos_vl_grl = pos_vl_grl,
    neg_tr_grl = neg_tr_grl,
    neg_vl_grl = neg_vl_grl
  )

  return(grl_list)
}


#' Keep Equal Number of Ranges
#'
#' This function randomly selects a specified number of genomic ranges from a
#' `GRanges` object to ensure equal representation.
#'
#' @param gr A `GRanges` object containing genomic ranges.
#' @param threshold An integer specifying the number of ranges to select.
#' @return A `GRanges` object containing the randomly selected ranges, sorted
#'   by genomic coordinates.
#' @importFrom GenomicRanges GRanges sort
#' @examples
#' gr <- GenomicRanges::GRanges(seqnames = "chr1", ranges = IRanges::IRanges(1:10, width = 1))
#' reduced_gr <- keep_equal_num(gr, 5)
#' @export
keep_equal_num <- function(gr, threshold) {
  random_indices <- sample(length(gr), threshold)
  gr_random <- gr[random_indices]
  gr_random <- GenomicRanges::sort(gr_random)
  return(gr_random)
}



#' Convert thick.start, thick.end, and thick.width to a thick IRanges Column
#'
#' This function takes a `GRangesList` object, where each `GRanges` contains 
#' columns `thick.start`, `thick.end`, and `thick.width` as metadata, and 
#' converts these into a single `thick` column of class `IRanges`.
#'
#' @param grl A `GRangesList` object with `thick.start`, `thick.end`, and 
#'   `thick.width` columns in each `GRanges`.
#' @return A `GRangesList` object where each `GRanges` contains a single 
#'   metadata column `thick` of class `IRanges`.
#' @importFrom GenomicRanges GRangesList
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors mcols
#' @examples
#' # Assuming `grl` is a GRangesList object with appropriate columns
#' modified_grl <- convert_to_thick_iranges(grl)
#' @export
convert_to_thick_iranges <- function(grl) {
  modified_grl <- lapply(grl, function(gr) {
    # Create IRanges object for the thick column
    thick_ir <- IRanges::IRanges(
      start = mcols(gr)$thick.start,
      end = mcols(gr)$thick.end
    )
    # Assign the IRanges object to the thick column
    S4Vectors::mcols(gr)$thick <- thick_ir
    # Optionally, drop other columns if needed
    S4Vectors::mcols(gr) <- S4Vectors::mcols(gr)[, "thick", drop = FALSE]
    return(gr)
  })

  return(GenomicRanges::GRangesList(modified_grl))
}
