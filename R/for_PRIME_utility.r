##' Create the output directory if it doesn't exist
#'
#' @param output_dir A character string specifying the path
#' to the output directory.
#'
#' @export
create_output_dir <- function(output_dir) {
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    message(sprintf("📁 Output directory created: %s", output_dir))
  } else {
    message(sprintf("📁 Output directory already exists: %s", output_dir))
  }
}

#' Set up logging to console or file
#'
#' @param log Logical or NULL. If NULL, log to console; otherwise,
#' log to file.
#' @param outdir Character. Path to the output directory
#' where logs will be stored.
#'
#' @return A connection or file path to be used as the logging target.
#' @export
setup_log_target <- function(log, outdir) {
  if (is.null(log)) {
    stdout()
  } else {
    log_dir <- file.path(outdir, "PRIMEloci_log")
    if (!dir.exists(log_dir)) {
      dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
      message(sprintf("📁 Log directory created: %s", log_dir))
    } else {
      message(sprintf("📁 Log directory already exists: %s", log_dir))
    }
    log_file <- file.path(log_dir, "PRIMEloci.log")
    message(sprintf("📝 Log file will be saved to: %s", log_file))
    log_file
  }
}

#' Set up temporary directory under output directory
#'
#' Creates a `PRIMEloci_tmp` folder under the given output directory
#' if it does not already exist.
#'
#' @param output_dir A character string specifying the path
#' to the output directory.
#'
#' @return A character string with the full path to the temporary directory.
#' @export
setup_tmp_dir <- function(output_dir) {
  tmp_dir <- file.path(output_dir, "PRIMEloci_tmp")
  if (!dir.exists(tmp_dir)) {
    dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
    message(sprintf("📁 Temporary directory created: %s", tmp_dir))
  } else {
    message(sprintf("📁 Temporary directory already exists: %s", tmp_dir))
  }
  tmp_dir
}

#' Extract CAGE Transcription Start Sites (CTSSs) from BigWig Files
#'
#' This function extracts CTSSs from CAGE data stored in BigWig files,
#' as specified in the provided design matrix. It allows for optional
#' genomic filtering based on a provided `Seqinfo` object and can exclude
#' or include specific chromosomes.
#'
#' @param dir_cage_bw A character string specifying the directory containing
#'   the BigWig files for both plus and minus strands.
#' @param design_matrix A data frame that must contain columns named
#'   \code{"BigWigPlus"} and \code{"BigWigMinus"}, which specify the paths
#'   to the respective BigWig files.
#' @param genome_info An optional \code{Seqinfo} object containing genome
#'   information. If provided, it allows for the use of `drop_chr`
#'   and `keep_chr` parameters.
#' @param drop_chr An optional character vector of chromosome names to be
#'   excluded from the analysis. Requires \code{genome_info} to be specified.
#' @param keep_chr An optional character vector of chromosome names to be
#'   included in the analysis. Requires \code{genome_info} to be specified.
#' @return A \code{CTSS} object containing the quantified and normalized
#'   counts of transcription start sites, including pooled counts.
#' @examples
#' \dontrun{
#' ctss_result <- get_ctss_from_bw("path/to/bigwig/files",
#'                                 design_matrix,
#'                                 genome_info,
#'                                 drop_chr = c("chrM"))
#' }
#' @import rtracklayer
#' @import GenomicRanges
#' @import CAGEfightR
#' @export
get_ctss_from_bw <- function(dir_cage_bw,
                             design_matrix,
                             genome_info = NULL,
                             drop_chr = NULL,
                             keep_chr = NULL) {

  # Create BigWigFileList objects for plus and minus strands
  bwplus <- rtracklayer::BigWigFileList(file.path(dir_cage_bw,
                                                  design_matrix$BigWigPlus))
  bwminus <- rtracklayer::BigWigFileList(file.path(dir_cage_bw,
                                                   design_matrix$BigWigMinus))
  names(bwplus) <- names(bwminus) <- rownames(design_matrix)

  # Handle genome_info, drop_chr, and keep_chr
  gn <- genome_info
  if (!is.null(genome_info)) {
    if (!is.null(drop_chr)) {
      gn <- GenomeInfoDb::dropSeqlevels(gn, drop_chr)
    }
    if (!is.null(keep_chr)) {
      gn <- GenomeInfoDb::keepSeqlevels(gn, keep_chr)
    }
  } else if (!is.null(drop_chr) || !is.null(keep_chr)) {
    stop("genome_info must be provided when using drop_chr or keep_chr.")
  }

  # Quantify CTSSs and calculate pooled counts
  orig_ctss <- CAGEfightR::quantifyCTSSs(plusStrand = bwplus,
                                         minusStrand = bwminus,
                                         design = design_matrix,
                                         genome = gn)
  orig_ctss <- CAGEfightR::calcPooled(orig_ctss, inputAssay = "counts")

  return(orig_ctss)
}
