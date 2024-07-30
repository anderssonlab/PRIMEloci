#' Import DHS Data from ENCODE
#'
#' This function reads DHS data from specified file paths, assuming the files 
#' are tab-delimited with columns "seqnames", "start", "end", and "name". It 
#' converts the start positions from 0-based to 1-based indexing and returns a 
#' list of data frames, each representing one file.
#'
#' @param dhs_paths_with_names A file path or a vector of file paths to the DHS 
#'   data files.
#' @return A list of data frames, each containing the DHS data from one file.
#' @importFrom readr cols col_double col_character read_delim
#' @importFrom purrr map
#' @importFrom dplyr mutate
#' @examples
#' # Example usage:
#' # Single file path:
#' # dhs_path <- "path/to/file1.txt"
#' # result <- import_dhs_encode(dhs_path)
#' #
#' # Multiple file paths:
#' # dhs_paths <- c("path/to/file1.txt", "path/to/file2.txt")
#' # result <- import_dhs_encode(dhs_paths)
#' @export
import_dhs_encode <- function(dhs_paths_with_names) {
  colnames <- c("seqnames", "start", "end", "name")
  coltypes <- readr::cols(.default = readr::col_double(),
                          seqnames = readr::col_character(),
                          name = readr::col_character())

  dhs_data <- purrr::map(dhs_paths_with_names, function(r) {
    df <- readr::read_delim(file = r, delim = "\t",
                            col_names = colnames, col_types = coltypes)
    df <- dplyr::mutate(df, start = start + 1) # Convert to 1-based indexing
    return(df)
  })

  return(dhs_data)
}

#' Filter Ranges Overlapping Non-reliable Genomic Regions
#'
#' Filters out ranges from a GRanges object that overlap with a specified
#' blacklist of non-reliable genomic regions.
#'
#' @param ranges A GRanges object containing the ranges to be filtered.
#' @param blacklist A GRanges object containing the blacklist of non-reliable
#'   genomic regions.
#' @return A GRanges object with ranges that do not overlap the blacklist.
#' @importFrom GenomicRanges findOverlaps
#' @importFrom S4Vectors queryHits
#' @examples
#' # Example usage:
#' # ranges <- GRanges(seqnames=Rle(c("chr1", "chr2")),
#' #                   ranges=IRanges(start=c(1, 3), end=c(5, 7)))
#' # blacklist <- GRanges(seqnames=Rle(c("chr1")),
#' #                      ranges=IRanges(start=c(2), end=c(4)))
#' # filtered_ranges <- filter_blacklist(ranges, blacklist)
#' @export
filter_blacklist <- function(ranges, blacklist) {
  overlapping_blacklist <- S4Vectors::queryHits(GenomicRanges::findOverlaps(ranges, blacklist))
  message("There are ", length(overlapping_blacklist), " ranges overlapping the blacklist")
  if (length(overlapping_blacklist) > 0) {
    ranges <- ranges[-overlapping_blacklist]
  }
  return(ranges)
}
