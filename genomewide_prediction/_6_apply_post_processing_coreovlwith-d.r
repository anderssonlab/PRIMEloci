#!/usr/bin/env Rscript

writeLines("\n### Running coreovlwith-d postprocess ###")

writeLines("\n# Importing R libraries..")
suppressPackageStartupMessages({
  library(argparse)
  library(GenomicRanges)
  library(parallel)
  library(tools)
  library(stringr)
  library(PRIMEloci)
  library(dplyr)
  library(data.table)
  library(IRanges)
})

# Create the argparse object
parser <- ArgumentParser()

parser$add_argument("-i", "--input_file", required = TRUE,
                    help = "Path to the input BED file.")
parser$add_argument("-t", "--score_threshold", type = "double", default = 0.7,
                    help = "Score threshold for filtering GRanges.")
parser$add_argument("-o", "--output_dir", default = NULL,
                    help = "Directory to save the output files.")
parser$add_argument("-m", "--use_max_cores", action = "store_true",
                    help = "Flag to use the maximum number of available cores.")
parser$add_argument("-d", "--score_diff", type = "double", default = 0.1,
                    help = "Score difference threshold for merging.")


# Parse arguments
args <- parser$parse_args()

# Define parameters based on parsed arguments
score_threshold <- args$score_threshold
score_diff <- args$score_diff
bed_file <- args$input_file
use_max_cores <- args$use_max_cores
output_dir <- args$output_dir

core_width <- 151

# Set the number of cores
num_cores <- if (use_max_cores) { parallel::detectCores() } else { 1 }

# Load and prepare data
bed <- load_bed_file(bed_file)
gr <- create_granges_from_bed(bed)

# Filter GRanges by score threshold
filtered_gr <- gr[gr$score >= score_threshold]
# Get a list of unique chromosomes
chr_list <- unique(as.character(seqnames(filtered_gr)))



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

      # Merge overlapping cores
      thick_vals <- c(thick_vals, thick_x)
      max_scores <- c(max_scores, score_x)
      merged_cores <- c(merged_cores, x)

      # Exclude merged cores
      core_gr <- core_gr[-1]
    }

    core_gr <- IRanges::subsetByOverlaps(core_gr,
                                         merged_cores,
                                         invert = TRUE)

  }

  # Attach metadata
  mcols(merged_cores)$thick <- thick_vals
  mcols(merged_cores)$max_score <- max_scores

  return(merged_cores)
}

write_granges_to_bed <- function(gr, output_dir, input_basename) {
  output_bed <- file.path(output_dir, paste0(input_basename,
                                             "_coreovlwith-d.bed"))
  bed_df <- as.data.frame(gr)

  if (!"thick" %in% colnames(mcols(gr)) || !"max_score" %in% colnames(mcols(gr))) { # nolint: line_length_linter.
    stop("'thick' or 'max_score' column not found in the metadata.")
  }

  bed_df$thick <- as.numeric(start(gr$thick))
  bed_df$start <- as.integer(bed_df$start - 1)  # Convert start to 0-based for BED format

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
                     col.names = TRUE)
  cat("Reduced GRanges object saved to", output_bed, "\n")
}



collapsed_gr_list <- mclapply(chr_list, function(chr) {

  tryCatch({

    chr_gr <- filtered_gr[GenomicRanges::seqnames(filtered_gr) == chr]

    core_gr <- GenomicRanges::resize(chr_gr, width = core_width, fix = "center")
    core_gr$thick <- start(core_gr) + floor(core_width / 2)

    final_merged_gr <- selective_merge_cores(core_gr, score_diff)

    return(final_merged_gr)
  }, error = function(e) {
    message(paste("Error processing chromosome", chr, ":", e$message))
    return(NULL)
  })
}, mc.cores = num_cores)

collapsed_gr <- do.call(c, collapsed_gr_list)
collapsed_gr <- sort(collapsed_gr)

if (!is.null(output_dir) && output_dir != FALSE) {
  input_basename <- tools::file_path_sans_ext(basename(bed_file))
  input_basename <- input_basename %>%
    stringr::str_replace_all("all", as.character(score_threshold)) %>%
    stringr::str_replace_all("[^[:alnum:]]", "_")
  write_granges_to_bed(collapsed_gr, output_dir, input_basename)
}

writeLines("Done!+++")
