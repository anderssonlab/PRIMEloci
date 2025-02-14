# Load necessary libraries
library(GenomicRanges)
library(readr)


shap_df <- readRDS("/Users/natsudanav/Documents/data_PRIMEloci_dev/SHAP/K562_endres_v2/K562_N_endres_v2_shap_value_df.rds")
subtnorm_df <- readRDS("/Users/natsudanav/Documents/data_PRIMEloci_dev/SHAP/K562_endres_v2/K562_N_endres_v2_profile_subtnorm_df.rds")
plus_count_df <- readRDS("/Users/natsudanav/Documents/data_PRIMEloci_dev/SHAP/K562_endres_v2/K562_N_endres_v2_profile_count_plus_df.rds")
minus_count_df <- readRDS("/Users/natsudanav/Documents/data_PRIMEloci_dev/SHAP/K562_endres_v2/K562_N_endres_v2_profile_count_minus_df.rds")
index <- readRDS("/Users/natsudanav/Documents/data_PRIMEloci_dev/SHAP/K562_endres_v2/K562_N_endres_v2_index.rds")

outname <- "K562_N_endres_v2"

head(shap_df)
head(index)

rownames(shap_df) <- index$rownames
rownames(subtnorm_df) <- index$rownames
rownames(plus_count_df) <- index$rownames
rownames(minus_count_df) <- index$rownames


parts <- strsplit(index$rownames, "[:;-]")
genomic_df <- do.call(rbind, parts)
genomic_df <- as.data.frame(genomic_df)
colnames(genomic_df) <- c("seqnames", "start", "end", "strand")

# Convert start and end to numeric
genomic_df$start <- as.numeric(genomic_df$start)
genomic_df$end <- as.numeric(genomic_df$end)

# Convert to GRanges object
gr <- GRanges(
  seqnames = genomic_df$seqnames,
  ranges = IRanges(start = genomic_df$start, end = genomic_df$end),
  strand = genomic_df$strand
)

# Print the GRanges object
gr

ls <- list()
ls$shap_value <- shap_df
ls$profile_subtnorm <- subtnorm_df
ls$profile_count_plus <- plus_count_df
ls$profile_count_minus <- minus_count_df
ls$gr <- gr

saveRDS(ls, file = paste0(outname, ".rds"))
