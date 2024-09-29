library(GenomicRanges)
library(IRanges)
library(parallel)
library(PRIMEloci)

# Define parameters
score_threshold <- 0.7
num_cores <- detectCores() - 1  # Adjust based on your system

# Load and prepare data
bed_file <- load_bed_file("/Users/natsudanav/Desktop/zmk214workingspace/data/resources/K562-on-PRIMEloci-sep-model_pred_all_profiles_subtnorm_tcs_K562_์.bed")
gr <- create_granges_from_bed(bed_file)

# Filter GRanges by score threshold
filtered_gr <- gr[gr$score >= score_threshold]

# Extract core function
extract_core <- function(gr) {
  start_core <- start(gr) + floor((width(gr) - 151) / 2)
  end_core <- start_core + 150  # 151bp total
  GRanges(seqnames(gr), IRanges(start = start_core, end = end_core), strand = strand(gr))
}

# Process by chromosome with metadata
process_by_chr <- function(chr, filtered_gr) {
  tryCatch({
    # Subset GRanges by chromosome
    chr_gr <- filtered_gr[seqnames(filtered_gr) == chr]
    core_gr <- extract_core(chr_gr)
    overlaps <- findOverlaps(core_gr)
    collapsed_ranges <- reduce(chr_gr[unique(queryHits(overlaps))])

    # Metadata for each region
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
    
    # Initialize empty metadata columns
    thick_vals <- numeric(length(collapsed_ranges))
    max_scores <- numeric(length(collapsed_ranges))
    all_scores <- character(length(collapsed_ranges))
    
    # Loop through each range in collapsed_ranges to extract metadata
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

# Run in parallel across chromosomes
chr_list <- unique(as.character(seqnames(filtered_gr)))
collapsed_gr_list <- mclapply(chr_list, process_by_chr, filtered_gr = filtered_gr, mc.cores = num_cores)

# Combine all results
collapsed_gr <- do.call(c, collapsed_gr_list)
collapsed_gr <- sort(collapsed_gr)


saveRDS(collapsed_gr, file = "collapsed_gr_N.rds")


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

  # Rearrange the columns to match the original BED file format
  #bed_df <- bed_df[, colnames(bed_file)]

  # Sort the data frame by chrom and chromStart
  bed_df <- bed_df[order(bed_df$chrom, bed_df$chromStart), ]

  # Write to BED file
  output_bed <- file.path(output_dir, paste0(input_basename, "_reduced.bed"))
  data.table::fwrite(bed_df, file = output_bed, sep = "\t", quote = FALSE, col.names = TRUE)

  cat("Reduced GRanges object saved to", output_bed, "\n")
}

save_granges_to_bed(collapsed_gr, ".", "firsttest_N")