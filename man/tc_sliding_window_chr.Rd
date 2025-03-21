% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tc_sld_window.r
\name{tc_sliding_window_chr}
\alias{tc_sliding_window_chr}
\title{Create Sliding Windows for One Chromosome}
\usage{
tc_sliding_window_chr(gr_per_chr, slide_by = 20, expand_by = 200)
}
\arguments{
\item{gr_per_chr}{A `GRanges` object representing a chromosome,
containing ranges of extended tag clusters (TC)
that have been uniformly extended by the same number of base pairs (bp).}

\item{slide_by}{An integer value specifying the size of
the sliding window step in base pairs (default: 20).}

\item{expand_by}{An integer value specifying
how much to expand the windows at both ends (default: 200 bp).}
}
\value{
A `GRanges` object containing the resized sliding windows
for the given chromosome.
}
\description{
This function generates sliding windows for a single chromosome,
starting from a reduced `GRanges` object.
The windows are resized by a specified amount,
and each window is expanded by a defined width.
The input `GRanges` object should represent
an extended Tag Cluster (TC) where all positions
have been uniformly extended by the same number of base pairs (bp).
}
\examples{
# Create a GRanges object
gr_chr <- GRanges(seqnames = "chr1", ranges = IRanges(start = c(100, 200), end = c(150, 250))) # nolint: line_length_linter.
# Generate sliding windows for the chromosome
sliding_windows <- tc_sliding_window_chr(gr_chr, slide_by = 20, expand_by = 200) # nolint: line_length_linter.
}
