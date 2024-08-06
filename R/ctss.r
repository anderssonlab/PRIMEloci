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
#'   counts of transcription start sites, including pooled counts
#'   and TPM values.
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

  # Quantify CTSSs, calculate pooled counts and TPM
  orig_ctss <- CAGEfightR::quantifyCTSSs(plusStrand = bwplus,
                                         minusStrand = bwminus,
                                         design = design_matrix,
                                         genome = gn)
  orig_ctss <- CAGEfightR::calcPooled(orig_ctss, inputAssay = "counts")
  orig_ctss <- CAGEfightR::calcTPM(orig_ctss,
                                   inputAssay = "counts",
                                   outputAssay = "TPM",
                                   totalTags = NULL,
                                   outputColumn = "totalTags")
  return(orig_ctss)
}
