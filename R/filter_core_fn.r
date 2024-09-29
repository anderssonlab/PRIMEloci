#' Extract the core region of a GRanges object.
#' 
#' @import GenomicRanges
#' @param gr A GRanges object.
#' @return A GRanges object representing the core region of each range (151 bp).
extract_core <- function(gr) {
  start_core <- start(gr) + floor((width(gr) - 151) / 2)
  end_core <- start_core + 150  # 151bp total
  GRanges(seqnames(gr), IRanges(start = start_core, end = end_core), strand = strand(gr))
}

#' Get metadata for each GRanges object.
#' 
#' @import GenomicRanges
#' @import IRanges
#' @param gr_range A single GRanges range.
#' @param filtered_gr A GRanges object with associated scores.
#' @return A list containing metadata: thick position, maximum score, and all scores.
get_metadata <- function(gr_range, filtered_gr) {
  overlap_hits <- subjectHits(findOverlaps(gr_range, filtered_gr))
  if (length(overlap_hits) > 0) {
    overlapping_scores <- filtered_gr$score[overlap_hits]
    max_score_region <- filtered_gr[overlap_hits][which.max(overlapping_scores)]
    max_score <- max(overlapping_scores)
    all_scores <- paste(overlapping_scores, collapse = ";")
    center_pos <- start(max_score_region) + floor(width(max_score_region) / 2)
    return(list(thick = center_pos, max_score = max_score, all_scores = all_scores))
  } else {
    return(list(thick = NA, max_score = NA, all_scores = NA))
  }
}

#' Process each chromosome and extract core regions and metadata.
#' 
#' @import GenomicRanges
#' @import IRanges
#' @param chr The chromosome name to process.
#' @param filtered_gr A filtered GRanges object with scores.
#' @return A GRanges object with metadata (thick, max_score, all_scores).
process_by_chr <- function(chr, filtered_gr) {
  tryCatch({
    # Subset GRanges by chromosome
    chr_gr <- filtered_gr[seqnames(filtered_gr) == chr]
    core_gr <- extract_core(chr_gr)
    overlaps <- findOverlaps(core_gr)
    collapsed_ranges <- reduce(chr_gr[unique(queryHits(overlaps))])

    # Initialize metadata columns
    thick_vals <- numeric(length(collapsed_ranges))
    max_scores <- numeric(length(collapsed_ranges))
    all_scores <- character(length(collapsed_ranges))
    
    # Loop through each collapsed range to extract metadata
    for (i in seq_along(collapsed_ranges)) {
      metadata <- get_metadata(collapsed_ranges[i], filtered_gr)
      thick_vals[i] <- metadata$thick
      max_scores[i] <- metadata$max_score
      all_scores[i] <- metadata$all_scores
    }
    
    # Assign metadata to collapsed ranges
    collapsed_ranges$thick <- thick_vals
    collapsed_ranges$max_score <- max_scores
    collapsed_ranges$all_scores <- all_scores
    
    return(collapsed_ranges)
  }, error = function(e) {
    message(paste("Error processing chromosome", chr, ":", e$message))
    return(NULL)  # Return NULL on error
  })
}

#' Save the GRanges object as both an RDS file and a BED file.
#' 
#' @import GenomicRanges
#' @import data.table
#' @param gr A GRanges object.
#' @param output_dir Directory where the output files will be saved.
#' @param input_basename The base name for the output files.
save_granges_to_bed <- function(gr, output_dir, input_basename) {
  # Save the selected GRanges object to an RDS file
  output_rds <- file.path(output_dir, paste0(input_basename, "_reduced.rds"))
  saveRDS(gr, file = output_rds)
  cat("Reduced GRanges object saved to", output_rds, "\n")

  # Convert GRanges object to a data frame
  bed_df <- as.data.frame(gr)

  # Select and rearrange the columns
  bed_df <- bed_df[, c("seqnames", "start", "end", "width", "strand", colnames(S4Vectors::mcols(gr)))]

  # Adjust the 'start' column to be 0-based by subtracting 1
  bed_df$start <- bed_df$start - 1

  # Rename the columns to match BED file format
  data.table::setnames(bed_df, c("seqnames", "start", "end", "strand"), c("chrom", "chromStart", "chromEnd", "strand"))

  # Sort the data frame by chrom and chromStart
  bed_df <- bed_df[order(bed_df$chrom, bed_df$chromStart), ]

  # Write to BED file
  output_bed <- file.path(output_dir, paste0(input_basename, "_reduced.bed"))
  data.table::fwrite(bed_df, file = output_bed, sep = "\t", quote = FALSE, col.names = TRUE)

  cat("Reduced GRanges object saved to", output_bed, "\n")
}
