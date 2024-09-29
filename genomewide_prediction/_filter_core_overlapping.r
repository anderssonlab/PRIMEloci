library(GenomicRanges)
library(IRanges)
library(parallel)
library(PRIMEloci)

library(S4Vectors)
library(data.table)


#' Extract the core region of a GRanges object.
#'
#' @import GenomicRanges
#' @param gr A GRanges object.
#' @return A GRanges object representing the core region of each range (151 bp).
extract_core <- function(gr, ext_core = 75) {
  start_core <- start(gr) + floor((width(gr) - (ext_core * 2 + 1)) / 2)
  end_core <- start_core + (ext_core * 2)  # 151bp total
  GRanges(seqnames(gr), IRanges(start = start_core,
                                end = end_core),
                                strand = strand(gr))
}

#' Get metadata for each GRanges object.
#'
#' @import GenomicRanges
#' @import IRanges
#' @param gr_range A single GRanges range.
#' @param filtered_gr A GRanges object with associated scores.
#' @return A list containing metadata: thick position,
#' maximum score, and all scores.
get_metadata <- function(gr_range, filtered_gr) {
  overlap_hits <- subjectHits(findOverlaps(gr_range, filtered_gr))
  if (length(overlap_hits) > 0) {
    overlapping_scores <- filtered_gr$score[overlap_hits]
    max_score_region <- filtered_gr[overlap_hits][which.max(overlapping_scores)]
    max_score <- max(overlapping_scores)
    all_scores <- paste(overlapping_scores, collapse = ";")
    center_pos <- start(max_score_region) + floor(width(max_score_region) / 2)
    return(list(thick = center_pos,
                max_score = max_score,
                all_scores = all_scores))
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
  bed_df <- bed_df[, c("seqnames", "start", "end", "width", "strand", colnames(S4Vectors::mcols(gr)))] # nolint: line_length_linter.

  # Adjust the 'start' column to be 0-based by subtracting 1
  bed_df$start <- bed_df$start - 1

  # Rename the columns to match BED file format
  data.table::setnames(bed_df,
                       c("seqnames", "start", "end", "strand"),
                       c("chrom", "chromStart", "chromEnd", "strand"))

  # Sort the data frame by chrom and chromStart
  bed_df <- bed_df[order(bed_df$chrom, bed_df$chromStart), ]

  # Write to BED file
  output_bed <- file.path(output_dir, paste0(input_basename, "_reduced.bed"))
  data.table::fwrite(bed_df,
                     file = output_bed,
                     sep = "\t",
                     quote = FALSE,
                     col.names = TRUE)

  cat("Reduced GRanges object saved to", output_bed, "\n")
}

# Load necessary libraries
# Assuming these are required packages for your script to run
library(GenomicRanges)
library(parallel)

#' Set the number of cores to use for parallel processing.
#'
#' This function determines the number of CPU cores
#' to use for parallel processing.
#' If the user specifies `use_max_cores = TRUE`,
#' the function will use the maximum
#' available cores minus one. Otherwise,
#' the user can specify the number of cores.
#'
#' @import parallel
#' @param use_max_cores Logical.
#' Whether to use the maximum number of cores minus one. Defaults to TRUE.
#' @param num_cores Integer. Number of cores specified by the user.
#' Ignored if `use_max_cores = TRUE`. Defaults to NULL.
#' @return Integer. The number of cores to use for parallel processing.
set_num_cores <- function(use_max_cores = TRUE, num_cores = NULL) {
  if (use_max_cores) {
    # Default to maximum cores - 1
    num_cores <- parallel::detectCores() - 1
    cat("Using", num_cores, "cores (maximum available minus one).\n")
  } else if (!is.null(num_cores)) {
    cat("Using", num_cores, "cores as specified.\n")
  } else {
    stop("Please specify the number of cores when use_max_cores = FALSE.")
  }
  return(num_cores)
}

# Define parameters
score_threshold <- 0.7

# Set the number of cores (use use_max_cores = TRUE by default)
num_cores <- set_num_cores(use_max_cores = TRUE)

# Load and prepare data
bed_file <- load_bed_file("/Users/natsudanav/Desktop/zmk214workingspace/data/resources/K562-on-PRIMEloci-sep-model_pred_all_profiles_subtnorm_tcs_K562_C.bed")
gr <- create_granges_from_bed(bed_file)

# Filter GRanges by score threshold
filtered_gr <- gr[gr$score >= score_threshold]

# Get a list of unique chromosomes
chr_list <- unique(as.character(seqnames(filtered_gr)))

# Run in parallel across chromosomes using the specified number of cores
collapsed_gr_list <- mclapply(chr_list,
                              process_by_chr,
                              filtered_gr = filtered_gr,
                              mc.cores = num_cores)

# Combine all results
collapsed_gr <- do.call(c, collapsed_gr_list)
collapsed_gr <- sort(collapsed_gr)

# Save the GRanges object to a file
save_granges_to_bed(collapsed_gr, ".", "firsttest_N")
saveRDS(collapsed_gr, file = "collapsed_gr_N.rds")

# If output_dir is specified, save the GRanges object to the directory
if (!is.null(output_dir) && output_dir != FALSE) {
  input_basename <- tools::file_path_sans_ext(basename(input_bed))
  save_granges_to_bed(selected_gr, output_dir, input_basename, bed_file)
}

writeLines("Done!")
