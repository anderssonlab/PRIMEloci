library(GenomicRanges)
library(PRIMEloci)
library(tidyverse)
library(ggplot2)

dist <- 200
ucsc_genome <- "hg38"
species <- "Homo sapiens"
genome_info <- GenomeInfoDb::keepStandardChromosomes(rtracklayer::SeqinfoForUCSCGenome(ucsc_genome), # nolint: line_length_linter.
                                                     species = species)
genome_info <- GenomeInfoDb::dropSeqlevels(genome_info, c("chrM", "chrY"))
standard_chromosomes <- GenomeInfoDb::seqnames(genome_info)

dhs_paths_with_names <- c("/Users/natsudanav/Desktop/zmk214workingspace/data/resources/E123-DNase.macs2.hg38.narrowPeak") # nolint: line_length_linter.
readin_dhs_data <- import_dhs_encode(dhs_paths_with_names)

dhs_ext <- lapply(readin_dhs_data, function(df) {
  df %>%
    extend_from_center_thick_df(.,
                                dist,
                                keep_same_length = TRUE,
                                group_exact_positions = TRUE) %>%
    dplyr::filter(seqnames %in% standard_chromosomes)
})
dhs_ext <- dhs_ext[[1]]
#dhs_thick <- calculate_thick_df(readin_dhs_data[[1]]) %>%
#  dplyr::filter(seqnames %in% standard_chromosomes) %>%
#  dplyr::arrange(seqnames, thick)

dhs_ext_gr <- GRanges(
  seqnames = dhs_ext$seqnames,
  ranges = IRanges(start = dhs_ext$start, end = dhs_ext$end),
  strand = "*"
)

dhs_thick_gr <- GRanges(
  seqnames = dhs_ext$seqnames,
  ranges = IRanges(start = dhs_ext$thick, end = dhs_ext$thick),
  strand = "*"
)


tc_grl <- readRDS("/Users/natsudanav/Desktop/zmk214workingspace/example/results/tc_grl.rds") # nolint: line_length_linter.

tc_gr <- tc_grl$K562_C3
tc_gr <- GenomeInfoDb::dropSeqlevels(tc_gr,
                                     c("chrM", "chrY"),
                                     pruning.mode = "coarse")
tc_thick_gr <- GRanges(
  seqnames = seqnames(tc_gr),
  ranges = tc_gr$thick,
  strand = strand(tc_gr),
  score = tc_gr$score
)
names(tc_thick_gr) <- NULL

tc_ovldhs <- IRanges::subsetByOverlaps(tc_thick_gr, dhs_ext_gr)

tc_ovldhs_pos <- tc_ovldhs[strand(tc_ovldhs) == "+" & tc_ovldhs$score >= 2]
tc_ovldhs_neg <- tc_ovldhs[strand(tc_ovldhs) == "-" & tc_ovldhs$score >= 2]

# Custom function to get the nearest features and filter by distance
get_nearest_features <- function(query_gr,
                                 reference_gr,
                                 direction = c("upstream", "downstream")) {

  # Choose direction based on input argument
  if (direction == "upstream") {
    nearest <- GenomicRanges::precede(query_gr, reference_gr)
  } else if (direction == "downstream") {
    nearest <- GenomicRanges::follow(query_gr, reference_gr)
  } else {
    stop("Invalid direction. Choose either 'upstream' or 'downstream'.")
  }

  # Filter out NA (no nearest found)
  nearest_non_na <- !is.na(nearest)
  query_filtered <- query_gr[nearest_non_na, ]
  nearest_filtered <- nearest[nearest_non_na]
  reference_filtered <- reference_gr[nearest_filtered, ]

  # Calculate distances between filtered query and reference
  distances <- GenomicRanges::distance(query_filtered, reference_filtered)

  return(distances)
}

tc_pos_downstream <- get_nearest_features(tc_ovldhs_pos,
                                          dhs_thick_gr,
                                          direction = "downstream")
tc_pos_downstream_200 <- tc_pos_downstream[tc_pos_downstream < 200]

tc_pos_upstream <- get_nearest_features(tc_ovldhs_pos,
                                        dhs_thick_gr,
                                        direction = "upstream")
tc_pos_upstream_200 <- tc_pos_upstream[tc_pos_upstream < 200]

tc_neg_downstream <- get_nearest_features(tc_ovldhs_neg,
                                          dhs_thick_gr,
                                          direction = "downstream")
tc_neg_downstream_200 <- tc_neg_downstream[tc_neg_downstream < 200]

tc_neg_upstream <- get_nearest_features(tc_ovldhs_neg,
                                        dhs_thick_gr,
                                        direction = "upstream")
tc_neg_upstream_200 <- tc_neg_upstream[tc_neg_upstream < 200]




file_name <- "Upstream_Downstream_Distances.pdf"
set_break <- 20
set_binwidth <- 200 / set_break


data_hist <- c(tc_neg_upstream_200, tc_pos_downstream_200)

# Create the histogram data without plotting to extract top bins
hist_data <- hist(data_hist, breaks = set_break, plot = FALSE)

# Get the counts and bin midpoints
bin_counts <- hist_data$counts
bin_mids <- hist_data$mids

# Sort the bin counts and select the top 5
top_5_indices <- order(bin_counts, decreasing = TRUE)[1:5]
top_5_counts <- bin_counts[top_5_indices]
top_5_mids <- bin_mids[top_5_indices]

print(top_5_mids)
print(top_5_counts)

# PLOTTING

# Combine the upstream (multiplied by -1) and downstream data
data_combined <- c((tc_neg_upstream_200 * -1), tc_pos_downstream_200)

# Create a data frame for easier plotting
data_df <- data.frame(
  distances = data_combined,
  direction = ifelse(data_combined < 0,
                     "Upstream", "Downstream")  # Label upstream and downstream
)

# Create the plot
plot <- ggplot(data_df, aes(x = distances, fill = direction)) +
  geom_histogram(binwidth = set_binwidth,
                 color = "black",
                 position = "identity",
                 alpha = 0.7) +
  scale_fill_manual(values = c("Upstream" = "lightblue",
                               "Downstream" = "orange")) +
  labs(
    title = "Upstream and Downstream Distance Distribution",
    x = "Distance (bp)",
    y = "Count"
  ) +
  theme_minimal() +
  theme(legend.title = element_blank())
print(plot)
# Save the plot to a PDF with the specified size and resolution
ggsave(file_name, plot = plot,
       width = 6, height = 8,
       units = "in", dpi = 500)
