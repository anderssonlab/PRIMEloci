library(GenomicRanges)
library(readr)
library(ggplot2)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(CAGEfightR)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

K562_C <- readRDS("/Users/natsudanav/Documents/data_PRIMEloci_dev/SHAP/K562_endres_v2/K562_C_endres_v2_shap-subtnorm-plus-minus-gr.rds")
K562_N <- readRDS("/Users/natsudanav/Documents/data_PRIMEloci_dev/SHAP/K562_endres_v2/K562_N_endres_v2_shap-subtnorm-plus-minus-gr.rds")

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

K562_C$gr <- CAGEfightR::assignTxType(K562_C$gr, txdb)
K562_N$gr <- CAGEfightR::assignTxType(K562_N$gr, txdb)

# Function convert GRanges to DataFrame and set the rownames
granges_to_df <- function(gr_obj) {
  df <- as.data.frame(gr_obj$gr)
  rownames(df) <- paste0(df$seqnames, ":",
                         df$start, "-",
                         df$end, ";",
                         df$strand)
  gr_obj$df <- df
  return(gr_obj)
}

K562_C <- granges_to_df(K562_C)
K562_N <- granges_to_df(K562_N)

# hclust_subtnorm <- function(profile_subtnorm) {
#   dist_matrix <- stats::dist(profile_subtnorm, method = "euclidean")
#   hc_wardd2 <- stats::hclust(dist_matrix, method = "ward.D2")
#   hc_average <- stats::hclust(dist_matrix, method = "average")
#   return(list(hc_wardd2, hc_average))
# }
# 
# hclust_fastcluster_subtnorm <- function(profile_subtnorm) {
#   dist_matrix <- stats::dist(profile_subtnorm, method = "euclidean")
#   hc_wardd2 <- fastcluster::hclust(dist_matrix, method = "ward.D2")
#   hc_average <- fastcluster::hclust(dist_matrix, method = "average")
#   return(list(hc_wardd2, hc_average))
# }

# names(K562_N)
# head(K562_N$profile_subtnorm)

# hc_K562_N <- hclust_fastcluster_subtnorm(K562_N$profile_subtnorm)
# saveRDS(hc_K562_N, "K562_N_endres_v2_fasthclust_subtnorm.rds")

# names(K562_C)
# head(K562_C$profile_subtnorm)

# hc_K562_C <- hclust_fastcluster_subtnorm(K562_C$profile_subtnorm)
# saveRDS(hc_K562_C, "K562_C_endres_v2_fasthclust_subtnorm.rds")

hc_K562_N <- readRDS("/Users/natsudanav/Documents/data_PRIMEloci_dev/SHAP/K562_endres_v2/K562_N_endres_v2_fasthclust_subtnorm.rds")
hc_K562_C <- readRDS("/Users/natsudanav/Documents/data_PRIMEloci_dev/SHAP/K562_endres_v2/K562_C_endres_v2_fasthclust_subtnorm.rds")

hc_K562_N
# > hc_K562_N
# [[1]]
# 
# Call:
# fastcluster::hclust(d = dist_matrix, method = "ward.D2")
# 
# Cluster method   : ward.D2 
# Distance         : euclidean 
# Number of objects: 116199 
# 
# 
# [[2]]
# 
# Call:
# fastcluster::hclust(d = dist_matrix, method = "average")
# 
# Cluster method   : average 
# Distance         : euclidean 
# Number of objects: 116199 


hc_C_wardD2 <- hc_K562_C[[1]]
hc_N_wardD2 <- hc_K562_N[[1]]

graphics::plot(hc_C_wardD2, main = "K562_C Hierarchical Clustering Dendrogram (euclidean distance, ward.D2 cluster)", xlab = "", sub = "")
graphics::plot(hc_N_wardD2, main = "K562_N Hierarchical Clustering Dendrogram (euclidean distance, ward.D2 cluster)", xlab = "", sub = "")

# Function to compute average profiles
compute_data_ave <- function(df_list) {

  if (length(df_list$shap_value) == 0) {
    warning("Empty data detected in compute_data_ave()")
    return(NULL)  # Return NULL if empty
  }

  data_ave <- data.frame(X = seq(-200, 200, length.out = 401))
  data_ave$SHAP <- colMeans(df_list$shap_value)
  data_ave$subtnorm <- colMeans(df_list$profile_subtnorm)
  data_ave$pos <- colMeans(df_list$profile_count_plus)
  data_ave$neg <- colMeans(df_list$profile_count_minus)
  return(data_ave)
}

# Function to generate plots including annotations
generate_plots <- function(data_ave, sample_name, annotation, num_data, df_list) {

  # Original SHAP and Subtnorm plots
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

  # New: Improved Annotation Bar Plot
  annotation_df <- as.data.frame(table(df_list$df$txType))
  colnames(annotation_df) <- c("Annotation", "Count")

  # Reorder annotations based on count (descending order)
  annotation_df$Annotation <- factor(annotation_df$Annotation,
                                     levels = annotation_df$Annotation[order(-annotation_df$Count)]) # nolint: line_length_linter.

  n4 <- ggplot(annotation_df, aes(x = Annotation,
                                  y = Count,
                                  fill = Annotation)) +
    geom_bar(stat = "identity", width = 0.5) + 
    scale_fill_manual(values = c(
      "promoter" = "#54854C", "proximal" = "#F3ECCA", "fiveUTR" = "#80B1D3",
      "threeUTR" = "#BEBEBE", "CDS" = "#000000", "exon" = "#C4D6A4",
      "intron" = "#D16227", "antisense" = "red", "intergenic" = "#FCDA70"
    )) +
    theme_minimal() +
    labs(title = paste("Annotation Distribution -",
                       sample_name, "- Cluster",
                       annotation),
         x = "Annotation Type", y = "Count") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")


  # Save all plots in one PDF file
  filename <- paste0("SHAP_Subtnorm_Profile_",
                     sample_name, "_",
                     annotation, ".pdf")

  if (file.exists(filename)) {
    file.remove(filename)
  }

  tryCatch({
    pdf(filename, width = 8, height = 12)
    grid.arrange(n1, n2, n3, n4, ncol = 1)  # Now includes annotation plot (n4)
    dev.off()
    message(paste("PDF saved as", filename))
  }, error = function(e) {
    message("Error while saving PDF for ", annotation, ": ", e$message)
  })
}


# Function to filter data frames based on annotations
filter_data_by_cluster <- function(data_list, k) {
  ref_df <- data_list$clusters
  filtered_list <- list()

  for (clus in seq(1, k)) {
    selected_rows <- rownames(ref_df[ref_df$clusters == clus, , drop = FALSE])
    filtered_list[[clus]] <- lapply(data_list, function(df) {
      if (is.data.frame(df)) {
        df[selected_rows, , drop = FALSE]
      }
    })
  }

  return(filtered_list)
}

k <- 13

clusters <- stats::cutree(hc_C_wardD2, k = k)
K562_C$clusters <- data.frame(clusters, row.names = names(clusters))
K562_C_filtered_df <- filter_data_by_cluster(K562_C, k)
for (clus in seq(1, k)) {
  df_list <- K562_C_filtered_df[[clus]]
  data_ave <- compute_data_ave(df_list)
  generate_plots(data_ave, "K562_C", clus, nrow(df_list$clusters), df_list)
}

l <- 12
clusters <- stats::cutree(hc_N_wardD2, k = l)
K562_N$clusters <- data.frame(clusters, row.names = names(clusters))
K562_N_filtered_df <- filter_data_by_cluster(K562_N, l)
for (clus in seq(1, l)) {
  df_list <- K562_N_filtered_df[[clus]]
  data_ave <- compute_data_ave(df_list)
  generate_plots(data_ave, "K562_N", clus, nrow(df_list$clusters), df_list)
}
