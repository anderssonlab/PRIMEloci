library(GenomicRanges)
library(readr)
library(ggplot2)
library(dendextend)
#library(ggfortify)
library(dplyr)
library(ggplot2)
#library(patchwork)
#library(gridExtra)
library(CAGEfightR)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

K562_C <- readRDS("/maps/projects/ralab/people/zmk214/K562_C_endres_v2_shap-subtnorm-plus-minus-gr.rds")
K562_N <- readRDS("/maps/projects/ralab/people/zmk214/K562_N_endres_v2_shap-subtnorm-plus-minus-gr.rds")

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



hclust_subtnorm <- function(profile_subtnorm) {
  dist_matrix <- stats::dist(profile_subtnorm, method = "euclidean")
  hc_wardd2 <- stats::hclust(dist_matrix, method = "ward.D2")
  hc_average <- stats::hclust(dist_matrix, method = "average")
  return(list(hc_wardd2, hc_average))
}

hclust_fastcluster_subtnorm <- function(profile_subtnorm) {
  dist_matrix <- stats::dist(profile_subtnorm, method = "euclidean")
  hc_wardd2 <- fastcluster::hclust(dist_matrix, method = "ward.D2")
  hc_average <- fastcluster::hclust(dist_matrix, method = "average")
  return(list(hc_wardd2, hc_average))
}

names(K562_N)
head(K562_N$profile_subtnorm)

hc_K562_N <- hclust_fastcluster_subtnorm(K562_N$profile_subtnorm)
saveRDS(hc_K562_N, "K562_N_endres_v2_fasthclust_subtnorm.rds")

names(K562_C)
head(K562_C$profile_subtnorm)

hc_K562_C <- hclust_fastcluster_subtnorm(K562_C$profile_subtnorm)
saveRDS(hc_K562_C, "K562_C_endres_v2_fasthclust_subtnorm.rds")

#library(cluster)
#hc <- cluster::diana(data)  # Performs divisive hierarchical clustering
#graphics::plot(hc)
