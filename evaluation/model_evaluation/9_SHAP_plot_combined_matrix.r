library(GenomicRanges)
library(readr)
library(ggplot2)
library(dplyr)
library(ggplot2)
library(patchwork)
library(gridExtra)
library(CAGEfightR)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

K562_C <- readRDS("/Users/natsudanav/Documents/data_PRIMEloci_dev/SHAP/K562_endres_v2/K562_C_endres_v2_shap-subtnorm-plus-minus-gr.rds")
K562_N <- readRDS("/Users/natsudanav/Documents/data_PRIMEloci_dev/SHAP/K562_endres_v2/K562_N_endres_v2_shap-subtnorm-plus-minus-gr.rds")

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

K562_C$gr <- CAGEfightR::assignTxType(K562_C$gr, txdb)
K562_N$gr <- CAGEfightR::assignTxType(K562_N$gr, txdb)

# Function convert GRanges to DataFrame and set the rownames
granges_to_df <- function(gr_obj, txdb) {
  df <- as.data.frame(gr_obj$gr)
  rownames(df) <- paste0(df$seqnames, ":",
                         df$start, "-",
                         df$end, ";",
                         df$strand)
  gr_obj$df <- df
  return(gr_obj)
}

K562_C <- granges_to_df(K562_C, txdb)
K562_N <- granges_to_df(K562_N, txdb)

# Function to compute average profiles
compute_data_ave <- function(df_list) {
  data_ave <- data.frame(X = seq(-200, 200, length.out = 401))
  data_ave$SHAP <- colMeans(df_list$shap_value)
  data_ave$subtnorm <- colMeans(df_list$profile_subtnorm)
  data_ave$pos <- colMeans(df_list$profile_count_plus)
  data_ave$neg <- colMeans(df_list$profile_count_minus)
  return(data_ave)
}

data_ave_C <- compute_data_ave(K562_C)
data_ave_N <- compute_data_ave(K562_N)


# Function to generate plots
generate_plots <- function(data_ave, sample_name, annotation, num_data) {

  n1 <- ggplot(data_ave, aes(x = X, y = SHAP)) +
    geom_bar(stat = "identity", fill = "#7771AF", alpha = 0.7) +
    labs(title = paste("Average SHAP Values (", sample_name, ") -",
                       annotation, "-", as.character(num_data), " rows"),
         x = "Position (bp)", y = "SHAP Value") +
    theme_minimal()

  n2 <- ggplot(data_ave, aes(x = X, y = subtnorm)) +
    geom_bar(stat = "identity", fill = "#3E8E74", alpha = 0.7) +
    labs(title = paste("Average Subtnorm Values (", sample_name, ") -",
                       annotation, "-", as.character(num_data), " rows"),
         x = "Position (bp)", y = "Subtnorm Value") +
    theme_minimal()

  n3 <- ggplot(data_ave, aes(x = X)) +
    geom_bar(aes(y = pos), stat = "identity", fill = "#264691", alpha = 0.7) +
    geom_bar(aes(y = -neg), stat = "identity", fill = "#DD553C", alpha = 0.7) +
    labs(title = paste("Average Profile Counts (", sample_name, ") - ",
                       annotation, " - ", as.character(num_data), " rows"),
         x = "Position (bp)", y = "Counts") +
    theme_minimal()

  # Ensure all plots are valid before saving
  if (!inherits(n1, "ggplot") || !inherits(n2, "ggplot") || !inherits(n3, "ggplot")) {
    message("One or more plots failed to generate. Skipping PDF saving for:",
            annotation)
    return(NULL)
  }

  # Save plots as PDF
  filename <- paste0("SHAP_Subtnorm_Profile_",
                     sample_name, "_",
                     annotation, ".pdf")

  if (file.exists(filename)) {
    file.remove(filename)
  }

  tryCatch({
    pdf(filename, width = 8, height = 12)
    grid.arrange(n1, n2, n3, ncol = 1)
    dev.off()
    message(paste("PDF saved as", filename))
  }, error = function(e) {
    message("Error while saving PDF for ", annotation, ": ", e$message)
  })

}

generate_plots(data_ave_C, "K562_C", "all", nrow(K562_C$df))
generate_plots(data_ave_N, "K562_N", "all", nrow(K562_N$df))

# Function to filter data frames based on annotations
filter_data_by_annotation <- function(data_list) {
  ref_df <- data_list$df
  annotation_categories <- unique(ref_df$txType)
  filtered_list <- list()

  for (ann in annotation_categories) {
    selected_rows <- rownames(ref_df[ref_df$txType == ann, ])
    filtered_list[[ann]] <- lapply(data_list, function(df) {
      if (is.data.frame(df)) {
        df[selected_rows, , drop = FALSE]
      }
    })
  }

  return(filtered_list)
}

K562_C_filtered_df <- filter_data_by_annotation(K562_C)
K562_N_filtered_df <- filter_data_by_annotation(K562_N)

for (ann in names(K562_C_filtered_df)) {
  df_list <- K562_C_filtered_df[[ann]]
  data_ave <- compute_data_ave(df_list)
  generate_plots(data_ave, "K562_C", ann, nrow(df_list$df))
}

for (ann in names(K562_N_filtered_df)) {
  df_list <- K562_N_filtered_df[[ann]]
  data_ave <- compute_data_ave(df_list)
  generate_plots(data_ave, "K562_N", ann, nrow(df_list$df))
}

