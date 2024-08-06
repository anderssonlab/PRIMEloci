#' Get Highest Non-Overlapping Ranges for a Chromosome
#'
#' This function takes a GenomicRanges::GRanges object
#' and returns the highest non-overlapping ranges
#' based on the scores in the metadata column.
#'
#' @param gr A GenomicRanges::GRanges object with a 'score' metadata column.
#' @return A GenomicRanges::GRanges object with
#' the highest non-overlapping ranges.
#' @import GenomicRanges
#' @import IRanges
#' @export
get_highest_non_overlap_chr <- function(gr) {

  gr <- GenomicRanges::sort(gr, by = ~score, decreasing = TRUE)
  result <- GenomicRanges::GRanges()

  while (length(gr) > 0) {
    highest <- gr[1]
    result <- c(result, highest)
    gr <- gr[!IRanges::overlapsAny(gr, highest)]
  }

  return(result)

}

#' Get Highest Non-Overlapping Ranges for All Chromosomes in Parallel
#'
#' This function takes a GenomicRanges::GRanges object
#' and returns the highest non-overlapping ranges
#' for all chromosomes, processed in parallel.
#'
#' @param gr A GenomicRanges::GRanges object with a 'score' metadata column.
#' @param num_cores The number of cores to use for parallel processing.
#' Default is parallel::detectCores() - 1.
#' @return A GenomicRanges::GRanges object with
#' the highest non-overlapping ranges for all chromosomes.
#'
#' @import GenomicRanges
#' @import foreach
#' @import doParallel
#' @import parallel
#' @export
get_highest_non_overlap <- function(gr,
                                    num_cores = parallel::detectCores() - 1) {

  gr_list <- GenomicRanges::split(gr, GenomicRanges::seqnames(gr))

  doParallel::registerDoParallel(cores = num_cores)

  results <- foreach::foreach(i = seq_along(gr_list), .combine = c) %dopar% {
    get_highest_non_overlap_chr(gr_list[[i]])
  }

  doParallel::stopImplicitCluster()

  results <- GenomicRanges::sort(results, by = ~start)

  return(results)

}
