# Load necessary libraries
library(GenomicRanges)
library(rtracklayer)
library(ggplot2)
library(reshape2)

# Define file paths
cage_file <- "results/all_te_prediction.bed"  # Path to your CAGE file
dhs_file <- "results/E116-DNase.macs2.hg38.narrowPeak"  # Path to your DHS file

# Load the CAGE data
cage_df <- read.table(cage_file, header=FALSE, sep="\t")
colnames(cage_df) <- c("chrom", "chromStart", "chromEnd", "name", "score", "strand", "sum_count", "annotation", "Class", "Dataset", "Augmentation", "CAGE", "Cell")

# Convert CAGE data to GRanges object
cage <- GRanges(seqnames=cage_df$chrom,
                ranges=IRanges(start=cage_df$chromStart, end=cage_df$chromEnd),
                sum_count=cage_df$sum_count)

# Load the DHS data (4-column file)
dhs_df <- read.table(dhs_file, header=FALSE, sep="\t")
colnames(dhs_df) <- c("chrom", "chromStart", "chromEnd", "name")

# Convert DHS data to GRanges object
dhs <- GRanges(seqnames=dhs_df$chrom,
               ranges=IRanges(start=dhs_df$chromStart, end=dhs_df$chromEnd))

# Perform overlaps between CAGE and DHS
overlaps <- findOverlaps(cage, dhs)

# Initialize DHS expression for all CAGE regions
cage$DHS_expression <- NA

# Aggregate `sum_count` for overlapping CAGE regions grouped by DHS regions
aggregated_dhs_expression <- tapply(cage$sum_count[queryHits(overlaps)],
                                     subjectHits(overlaps), sum, na.rm=TRUE)

# Assign the aggregated DHS expression back to the respective DHS regions
cage$DHS_expression[queryHits(overlaps)] <- aggregated_dhs_expression[subjectHits(overlaps)]

# Bin DHS expression into 10 levels
cage$dhs_expression_bin <- cut(cage$DHS_expression, breaks=10, labels=FALSE)

# Compute recall for each bin
recall_summary <- aggregate(cage$DHS_expression, 
                             by=list(bin=cage$dhs_expression_bin), 
                             FUN=function(x) length(which(!is.na(x))) / length(cage))
colnames(recall_summary) <- c("DHS_expression", "Recall")

# Plot Recall vs. DHS Expression
ggplot(recall_summary, aes(x=DHS_expression, y=Recall)) +
  geom_line() +
  labs(title="Recall vs. DHS Expression", x="DHS Expression", y="Recall")

# Prepare data for heatmap
heatmap_df <- data.frame(CAGE_sum_count=cage$sum_count,
                         DHS_expression=cage$DHS_expression)
heatmap_data <- melt(table(cut(heatmap_df$CAGE_sum_count, breaks=10),
                           cut(heatmap_df$DHS_expression, breaks=10)))

# Plot Heatmap of Sum Count vs. DHS Expression
ggplot(heatmap_data, aes(x=Var1, y=Var2, fill=value)) +
  geom_tile() +
  scale_fill_gradient(low="blue", high="red") +
  labs(title="Heatmap: Sum Count vs. DHS Expression", x="CAGE Sum Count", y="DHS Expression", fill="Count")
