#' Create Sliding Windows for One Chromosome
#'
#' This function generates sliding windows for a single chromosome,
#' starting from a reduced `GRanges` object.
#' The windows are resized by a specified amount,
#' and each window is expanded by a defined width.
#' The input `GRanges` object should represent
#' an extended Tag Cluster (TC) where all positions
#' have been uniformly extended by the same number of base pairs (bp).
#'
#' @param gr_per_chr A `GRanges` object representing a chromosome,
#' containing ranges of extended tag clusters (TC)
#' that have been uniformly extended by the same number of base pairs (bp).
#' @param slide_by An integer value specifying the size of
#' the sliding window step in base pairs (default: 20).
#' @param expand_by An integer value specifying
#' how much to expand the windows at both ends (default: 200 bp).
#' @return A `GRanges` object containing the resized sliding windows
#' for the given chromosome.
#' @examples
#' # Create a GRanges object
#' gr_chr <- GRanges(seqnames = "chr1", ranges = IRanges(start = c(100, 200), end = c(150, 250))) # nolint: line_length_linter.
#' # Generate sliding windows for the chromosome
#' sliding_windows <- tc_sliding_window_chr(gr_chr, slide_by = 20, expand_by = 200) # nolint: line_length_linter.
#' @import GenomicRanges
#' @import IRanges
tc_sliding_window_chr <- function(gr_per_chr, slide_by = 20, expand_by = 200) { # nolint: line_length_linter.

  # 1) Reduce overlaps within this chromosome
  collapsed_granges <- GenomicRanges::reduce(gr_per_chr)

  # 2) Generate sliding windows for all ranges in the collapsed GRanges
  sliding_granges_list <- base::lapply(seq_along(collapsed_granges), function(i) { # nolint: line_length_linter.
    start_pos <- GenomicRanges::start(collapsed_granges[i])
    end_pos <- GenomicRanges::end(collapsed_granges[i])

    # Adjust the sliding window start and end
    adjusted_start <- start_pos + expand_by - slide_by
    adjusted_end <- end_pos - expand_by + slide_by

    # Create a sequence of positions for sliding windows
    sliding_positions <- base::seq(adjusted_start, adjusted_end, by = slide_by)

    # Create a GRanges object for the sliding windows (width = 1)
    sliding_granges <- GenomicRanges::GRanges(
      seqnames = GenomicRanges::seqnames(collapsed_granges[i]),
      ranges = IRanges::IRanges(start = sliding_positions, width = 1),
      strand = GenomicRanges::strand(collapsed_granges[i])
    )

    # Store the original ranges (before expansion) as IRanges
    # in a "thick" metadata column
    sliding_granges$thick <- IRanges::IRanges(start = sliding_positions,
                                              width = 1)

    # Expand sliding_granges directly to save memory
    sliding_granges <- GenomicRanges::resize(sliding_granges,
                                             width = GenomicRanges::width(sliding_granges) + expand_by, # nolint: line_length_linter.
                                             fix = "start")
    sliding_granges <- GenomicRanges::resize(sliding_granges,
                                             width = GenomicRanges::width(sliding_granges) + expand_by, # nolint: line_length_linter.
                                             fix = "end")

    return(sliding_granges)
  })

  # 3) Combine the list of GRanges objects into a single GRanges object
  sliding_granges <- base::do.call(c, sliding_granges_list)

  # Return the sliding GRanges object
  return(sliding_granges)
}



#' Create Sliding Windows Genome-Wide from Tag Clusters and Combine the Results
#'
#' This function applies sliding window generation
#' for each chromosome in a `GRanges` object
#' representing extended tag clusters (TCs)
#' and processes all chromosomes in parallel.
#' The input `GRanges` object must represent tag clusters
#' that have been uniformly extended by the same number of base pairs.
#' The sliding windows are combined into a `GRangesList`,
#' which is then flattened into a single `GRanges` object.
#'
#' @param granges_obj A `GRanges` object representing
#' the genome-wide extended tag clusters (TCs) where each
#' range has been uniformly extended by the same number of base pairs (bp).
#' @param slide_by An integer value specifying the size of
#' the sliding window step in base pairs (default: 20).
#' @param expand_by An integer value specifying
#' how much to expand the windows at both ends (default: 200 bp).
#' @param use_max_cores Logical, if `TRUE` will use almost
#' all available CPU cores for parallel processing (default: TRUE).
#' @return A `GRanges` object containing the combined sliding windows
#' for all chromosomes in the input.
#' @examples
#' # Create a GRanges object
#' gr <- GRanges(seqnames = c("chr1", "chr2"), ranges = IRanges(start = c(100, 500), end = c(150, 600))) # nolint: line_length_linter.
#' # Generate sliding windows for the entire genome
#' tc_sliding_windows <- tc_sliding_window(gr, slide_by = 20, expand_by = 200) # nolint: line_length_linter.
#' @import GenomicRanges
#' @import IRanges
#' @import parallel
#' @export
tc_sliding_window <- function(granges_obj,
                              slide_by = 20,
                              expand_by = 200,
                              use_max_cores = TRUE) {

  # 1) Ignore strand by setting all strands to "*"
  GenomicRanges::strand(granges_obj) <- "*"

  # 2) Split the GRanges object by chromosome (seqnames)
  gr_by_chr <- base::split(granges_obj, GenomicRanges::seqnames(granges_obj))

  # 3) Determine the number of cores to use
  if (use_max_cores) {
    # Use almost full capacity, leaving 1 core free
    mc_cores <- parallel::detectCores() - 1
  } else {
    # Default to 2 cores if not using full capacity
    mc_cores <- 2
  }

  # 4) Parallelize the processing over each chromosome using mclapply
  result_list <- parallel::mclapply(gr_by_chr,
                                    tc_sliding_window_chr,
                                    slide_by = slide_by,
                                    expand_by = expand_by,
                                    mc.cores = mc_cores)

  # 5) Convert the result into a GRangesList
  result_grl <- GenomicRanges::GRangesList(result_list)

  # 6) Unlist the GRangesList to create a single GRanges object
  result_gr <- unlist(result_grl)

  # Return the final unlisted GRanges object
  return(result_gr)
}
