#!/usr/bin/env Rscript
# PCA plot for the clusters indentified by Souporcell

# Usage: 
# Rscript pca_souporcell_clusters <souporcell_clusters.tsv> <output_file>

options(stringsAsFactor = FALSE)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  stop("Not enough arguments, requires two!\nUsage:\n  Rscript plot_pca_v2.r souporcell_clusters.tsv output.pdf")
}

library(tidyverse)
library(magrittr)

# Load the cluster table
df <- readr::read_table(args[1])

# If the cell is assigned to multiple cluster, set the assignment to "Undetermined"
cluster_tab <- dplyr::mutate(
  df, final_assignment = dplyr::if_else(stringr::str_detect(assignment, "/"), "Undetermined", assignment)
) 

# Obtain an vector containing assignments and named by cell barcode
assignment_tab <- dplyr::pull(cluster_tab, final_assignment, barcode) 

# Do PCA for all available clusters
pca <- tibble::column_to_rownames(cluster_tab, "barcode") %>%
  dplyr::select(dplyr::starts_with("cluster")) %>%
  (function(tab) { tab / rowMeans(tab) }) %>%
  prcomp(cneter = TRUE, scale. = TRUE) %>%
  (function(tab) { as.data.frame(tab$x) %>% dplyr::mutate(colors = assignment_tab[rownames(.)]) })

# A PCA plot and save it to the disk
pca_plot <- ggplot() +
  geom_point(aes(x = PC1, y = PC2), dplyr::filter(pca, colors == "Undetermined"), color = "gray", alpha = 0.5, size = 0.5) +
  geom_point(aes(x = PC1, y = PC2, color = colors), dplyr::filter(pca, colors != "Undetermined"), alpha = 0.5, size = 0.5)  +
  theme_bw() +
  theme(panel.grid = element_blank(), aspect.ratio = 1) +
  labs(x = "PC1", y = "PC2", color = "Cluster")
ggsave(args[2], plot = pca_plot, width = 5, height = 5)
