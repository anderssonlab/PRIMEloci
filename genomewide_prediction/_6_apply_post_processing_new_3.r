#!/usr/bin/env Rscript

writeLines("\n### Running NEW POSTPROCESS ###")

writeLines("\n# Importing R libraries..")
suppressPackageStartupMessages({
  library(argparse)
  library(GenomicRanges)
  library(parallel)
  library(tools)
  library(stringr)
  library(PRIMEloci)
})

library(GenomicRanges)
library(IRanges)
library(S4Vectors)

# Create the argparse object
parser <- ArgumentParser(description = "Process GRanges with selective core merging based on score and width thresholds.") # nolint: line_length_linter.
parser$add_argument("-i", "--input_file", required = TRUE,
                    help = "Path to the input BED file.")
parser$add_argument("-t", "--score_threshold", type = "double", default = 0.7,
                    help = "Score threshold for filtering GRanges.")
parser$add_argument("-d", "--score_diff", type = "double", default = 0.1,
                    help = "Score difference threshold for merging.")
parser$add_argument("-w", "--max_width", type = "integer", default = 1000,
                    help = "Maximum width for merged regions.")
parser$add_argument("-m", "--use_max_cores", action = "store_true",
                    default = TRUE,
                    help = "Flag to use the maximum number of available cores.")
parser$add_argument("-o", "--output_dir", default = NULL,
                    help = "Directory to save the output files.")

# Parse arguments
args <- parser$parse_args()

# Define parameters based on parsed arguments
score_threshold <- as.numeric(args$score_threshold)
score_diff <- as.numeric(args$score_diff)
max_width <- as.numeric(args$max_width)
bed_file <- args$input_file
use_max_cores <- args$use_max_cores
output_dir <- args$output_dir

# Set the number of cores
num_cores <- if (use_max_cores) { parallel::detectCores() } else { 1 }

# Load and prepare data
bed <- load_bed_file(bed_file)
gr <- create_granges_from_bed(bed)

# Filter GRanges by score threshold
filtered_gr <- gr[gr$score >= score_threshold]
print(head(filtered_gr))
# Get a list of unique chromosomes
chr_list <- unique(as.character(seqnames(filtered_gr)))
print(chr_list)

selective_merge_cores <- function(core_gr, score_diff, maxgap = 150, max_width) {
  # Sort by score (highest to lowest)
  core_gr <- core_gr[order(-core_gr$score)]
  
  # Initialize empty GRanges for final results
  merged_cores <- GenomicRanges::GRanges()
  
  # Track used regions to avoid reprocessing
  used <- logical(length(core_gr))
  
  for (i in seq_along(core_gr)) {
    # Skip if region is already part of a merge
    if (used[i]) next
    
    # Select the current region
    current_region <- core_gr[i]
    
    # Find regions within maxgap
    nearby_indices <- which(start(core_gr) >= start(current_region) - maxgap &
                            start(core_gr) <= end(current_region) + maxgap &
                            !used)
    
    # Filter nearby regions based on score_diff
    merge_candidates <- core_gr[nearby_indices][core_gr[nearby_indices]$score >= current_region$score - score_diff]
    
    # If there are multiple merge candidates, merge them
    if (length(merge_candidates) > 1) {
      # Merge candidates manually
      start_pos <- min(start(merge_candidates))
      end_pos <- max(end(merge_candidates))
      
      # Create the merged region without width constraints
      merged_region <- GRanges(seqnames = seqnames(current_region),
                               ranges = IRanges(start = start_pos, end = end_pos),
                               strand = strand(current_region),
                               score = max(mcols(merge_candidates)$score))  # Keep highest score

      # Add metadata for 'thick' as midpoint
      mcols(merged_region)$thick <- IRanges::IRanges(start = start_pos + floor((end_pos - start_pos) / 2), width = 1)
      
      # Append to final result
      merged_cores <- c(merged_cores, merged_region)
      
      # Mark candidates as used
      used[nearby_indices] <- TRUE
    } else {
      # No valid merge candidates; add the current region directly
      merged_cores <- c(merged_cores, current_region)
      used[i] <- TRUE
    }
  }
  
  return(merged_cores)
}



library(GenomicRanges)
library(IRanges)
library(data.table)

# Function to write a GRanges object to a BED file with input_basename in the file name
write_granges_to_bed <- function(gr, output_dir, input_basename) {
  # Construct the output file name based on input_basename
  output_bed <- file.path(output_dir, paste0(input_basename, "_new2.bed"))
  
  # Convert GRanges object to data frame
  bed_df <- as.data.frame(gr)

  # Ensure the 'thick' and 'max_score' columns exist in the metadata
  if (!"thick" %in% colnames(mcols(gr))) {
    stop("'thick' column not found in the metadata.")
  }
  if (!"max_score" %in% colnames(mcols(gr))) {
    stop("'max_score' column not found in the metadata.")
  }

  # Convert 'thick' from IRanges to numeric
  bed_df$thick <- as.numeric(start(gr$thick))

  # Select the columns for BED format (seqnames, start, end, strand, thickStart, thickEnd, max_score)
  bed_df <- bed_df[, c("seqnames", "start", "end", "strand", "thick", "max_score")]

  # Adjust the 'start' column to be 0-based for BED format
  bed_df$start <- bed_df$start - 1

  # Create thickStart and thickEnd columns using 'thick'
  bed_df$thickStart <- bed_df$thick - 1  # 0-based thickStart for BED format
  bed_df$thickEnd <- bed_df$thick        # thickEnd is the 1-based thick position

  # Remove the 'thick' column as it's been split into thickStart and thickEnd
  bed_df <- bed_df[, !colnames(bed_df) %in% "thick"]

  # Generate a name column based on row index, using 'corereduced_XX'
  bed_df$name <- paste0("ovlcorereduced_", seq_len(nrow(bed_df)))

  # Reorder the columns to match the expected order
  desired_cols <- c("seqnames", "start", "end", "name", "max_score", "strand", "thickStart", "thickEnd")
  bed_df <- bed_df[, desired_cols]

  # Rename the columns to match BED format
  data.table::setnames(bed_df,
                       c("seqnames", "start", "end", "strand"),
                       c("chrom", "chromStart", "chromEnd", "strand"))

  # Write to BED file
  output_bed <- file.path(output_dir,
                          paste0(input_basename,
                                 "_new3.bed"))
  data.table::fwrite(bed_df,
                     file = output_bed,
                     sep = "\t",
                     quote = FALSE,
                     col.names = FALSE)

  cat("Reduced GRanges object saved to", output_bed, "\n")
}





# Run in parallel across chromosomes using the specified number of cores
collapsed_gr_list <- mclapply(chr_list, function(chr) {
  tryCatch({
    # Subset GRanges by chromosome
    chr_gr <- filtered_gr[GenomicRanges::seqnames(filtered_gr) == chr]
    core_gr <- GenomicRanges::resize(chr_gr, width = 151, fix = "center")
    print(head(core_gr))
    # Apply selective merging on overlapping sets of cores
    final_merged_gr <- selective_merge_cores(core_gr, score_diff, max_width)
    print(head(final_merged_gr))

    return(final_merged_gr)

  }, error = function(e) {
    message(paste("Error processing chromosome", chr, ":", e$message))
    return(NULL)
  })
}, mc.cores = num_cores)

# Combine all results
collapsed_gr <- do.call(c, collapsed_gr_list)
collapsed_gr <- sort(collapsed_gr)

# If output_dir is specified, save the GRanges object to the directory
if (!is.null(output_dir) && output_dir != FALSE) {
  input_basename <- tools::file_path_sans_ext(basename(bed_file))
  input_basename <- input_basename %>%
    stringr::str_replace_all("all", as.character(score_threshold)) %>%
    stringr::str_replace_all("[^[:alnum:]]", "_")

  write_granges_to_bed(collapsed_gr, output_dir, input_basename)
}

writeLines("Done!")
