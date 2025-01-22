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
    name <- paste0("PRIMEloci-", seq_len(nrow(bed_df)))
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