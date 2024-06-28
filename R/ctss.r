#' Extract CTSSs from BigWig files
#'
#' This function processes CAGE data from BigWig files specified in the
#' design matrix, optionally considering genomic information and allowing
#' for chromosome level filtering.
#'
#' @param dir_cage_bw Directory containing the BigWig files.
#' @param design_matrix Data frame with BigWig file information. It should
#'   contain columns named \code{"BigWigPlus"} and \code{"BigWigMinus"}.
#' @param genome_info Optional \code{Seqinfo} object
#' containing genome information.
#' @param drop_chr Optional character vector of chromosomes to drop.
#' @param keep_chr Optional character vector of chromosomes to keep.
#'
#' @return A processed CTSS object with quantified and normalized counts.
#'
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
get_ctss_from_bw <- function(dir_cage_bw, design_matrix,
                             genome_info = NULL,
                             drop_chr = NULL, keep_chr = NULL) {

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
      gn <- GenomicRanges::dropSeqlevels(gn, drop_chr)
    }
    if (!is.null(keep_chr)) {
      gn <- GenomicRanges::keepSeqlevels(gn, keep_chr)
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
