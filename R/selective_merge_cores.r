#' Selectively merge overlapping cores based on score difference.
#'
#' This function takes a GRanges object of genomic cores, sorts them
#' by descending score, and merges overlapping cores if their score
#' difference is within a specified threshold.
#'
#' @param core_gr A GRanges object containing genomic core ranges.
#' @param score_diff A numeric value specifying the maximum allowed
#' score difference for merging overlapping cores.
#' @return A GRanges object with merged core regions and added
#' metadata, including the thick position and maximum score.
#' @importFrom GenomicRanges GRanges reduce findOverlaps
#' @importFrom IRanges IRanges subsetByOverlaps
#' @importFrom S4Vectors mcols subjectHits
#' @export
selective_merge_cores <- function(core_gr, score_diff) {
  # Sort cores by descending score
  core_gr <- core_gr[order(-core_gr$score)]

  # Initialize final GRanges and metadata containers
  merged_cores <- GenomicRanges::GRanges()
  thick_vals <- IRanges::IRanges()
  max_scores <- numeric(0)

  while (length(core_gr) > 0) {
    # Take the highest-ranked core
    x <- core_gr[1]
    score_x <- x$score
    thick_x <- x$thick

    # Find and filter overlapping cores
    overlaps <- GenomicRanges::findOverlaps(x, core_gr)
    overlap_set <- core_gr[S4Vectors::subjectHits(overlaps)]
    merge_candidates <- overlap_set[overlap_set$score >= score_x - score_diff]

    if (length(merge_candidates) > 1) {
      # Merge overlapping cores
      merged_region <- GenomicRanges::reduce(merge_candidates)
      thick_vals <- c(thick_vals, thick_x)
      max_scores <- c(max_scores, score_x)
      merged_cores <- c(merged_cores, merged_region)

      # Exclude merged cores
      core_gr <- core_gr[!(seqnames(core_gr) %in% seqnames(overlap_set) &
                             start(core_gr) %in% start(overlap_set) &
                             end(core_gr) %in% end(overlap_set))]
    } else {
      thick_vals <- c(thick_vals, thick_x)
      max_scores <- c(max_scores, score_x)
      merged_cores <- c(merged_cores, x)

      # Exclude merged cores
      core_gr <- core_gr[-1]
    }

    core_gr <- IRanges::subsetByOverlaps(core_gr, merged_cores, invert = TRUE)
  }

  # Attach metadata
  mcols(merged_cores)$thick <- thick_vals
  mcols(merged_cores)$max_score <- max_scores

  return(merged_cores)
}