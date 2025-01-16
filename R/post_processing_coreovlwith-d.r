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
#' @importFrom S4Vectors mcols
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
    overlap_set <- core_gr[subjectHits(overlaps)]
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


#' Write a GRanges object to a BED file for coreovlwith-d.
#'
#' This function converts a GRanges object to a data frame
#' and writes it to a BED file, including metadata such as
#' thick position and maximum score.
#'
#' @param gr A GRanges object containing genomic ranges and metadata.
#' @param output_dir The directory where the BED file will be saved.
#' @param input_basename The base name for the output BED file.
#' @param score_diff A numeric value representing the score difference
#' used for naming the BED file.
#' @return None. The function writes a BED file to the specified
#' output directory.
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors mcols
#' @importFrom data.table fwrite
#' @export
write_granges_to_bed_coreovlwithd <- function(gr,
                                              output_dir,
                                              input_basename,
                                              score_diff) {

  output_bed <- file.path(output_dir, paste0(input_basename,
                                             "coreovlwith-d",
                                             str(score_diff),
                                             ".bed"))
  bed_df <- as.data.frame(gr)

  if (!"thick" %in% colnames(mcols(gr)) || !"max_score" %in% colnames(mcols(gr))) { # nolint: line_length_linter.
    stop("'thick' or 'max_score' column not found in the metadata.")
  }

  bed_df$thick <- as.numeric(start(gr$thick))
  bed_df$start <- bed_df$start - 1  # Convert start to 0-based for BED format

  bed_df <- within(bed_df, {
    thickStart <- thick - 1
    thickEnd <- thick
    name <- paste0("coreovlwith-d", seq_len(nrow(bed_df)))
  })

  # Reorder and rename columns
  bed_df <- bed_df[, c("seqnames", "start", "end",
                       "name", "max_score", "strand",
                       "thickStart", "thickEnd")]
  data.table::setnames(bed_df,
                       c("seqnames", "start", "end", "strand"),
                       c("chrom", "chromStart", "chromEnd", "strand"))

  # Write to BED file
  data.table::fwrite(bed_df,
                     file = output_bed,
                     sep = "\t",
                     quote = FALSE,
                     col.names = FALSE)
  cat("Reduced GRanges object saved to", output_bed, "\n")
}