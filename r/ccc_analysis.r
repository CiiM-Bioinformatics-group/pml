#!/usr/bin/env Rscript
# File: deg_analysis.r
# Author: Zhenhua Zhang
# E-mail: zhenhua.zhang217@gmail.com
# Created: Jun 28, 2023
# Updated: Nov 26, 2023

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)

  library(Seurat)
  library(NMF)
  library(CellChat)
  library(patchwork)
})

set.seed(31415926)

rm(list = ls()); gc()
overwrite <- TRUE
proj_dir <- "~/Documents/projects/wp_pml"
plot_dir <- file.path(proj_dir, "outputs/analysis/cell_cell_communication/plots")
object_dir <- file.path(proj_dir, "outputs/analysis/cell_cell_communication/objects")

tar_cell_types <- c("Mono", "CD4 T", "CD8 T", "B", "NK", "DC")

per_omics <- "pml_citeseq"
# per_omics <- "pml_rnaseq"

response <- "Responder"; timepoint <- "BL"
# response <- "Responder"; timepoint <- "6W"
# response <- "Responder"; timepoint <- "3M"
# response <- "Non-responder"; timepoint <- "BL"
# response <- "Non-responder"; timepoint <- "3M"
save_token <- paste("ccc", per_omics, response, timepoint, sep = ".")


#
## Cell communications by CellChat
#
# CellChat DB of human
cc_db.use <- CellChatDB.human %>% subsetDB(search = "Secreted Signaling")

# Create CellChat object
pbmc_int <- file.path(proj_dir, "outputs/analysis/integrated/objects", paste0(per_omics, ".rds")) %>% readRDS()
tar_cells <- pbmc_int@meta.data %>%
  dplyr::filter(Response == response, Timepoint == timepoint, predicted.celltype.l1 %in% tar_cell_types) %>%
  rownames()

sub_pbmc <- pbmc_int[, tar_cells]
data_input <- GetAssayData(sub_pbmc, assay = "SCT", slot = "data") # normalized data matrix
meta_data <- sub_pbmc@meta.data %>%
  dplyr::select(PatientID = Vireo_assignment, group = predicted.celltype.l1, response = Response, Timepoint = Timepoint)

cc_obj <- createCellChat(object = data_input, meta = meta_data, group.by = "group")
cc_obj <- setIdent(cc_obj, ident.use = "group")

# CellChat analysis
cc_obj@DB <- cc_db.use
cc_obj <- subsetData(cc_obj)
cc_obj <- identifyOverExpressedGenes(cc_obj)
cc_obj <- identifyOverExpressedInteractions(cc_obj)
cc_obj <- computeCommunProb(cc_obj, population.size = TRUE)
cc_obj <- filterCommunication(cc_obj, min.cells = 10)
cc_obj <- computeCommunProbPathway(cc_obj)
cc_obj <- aggregateNet(cc_obj)
cc_obj <- netAnalysis_computeCentrality(cc_obj, slot.name = "netP")


# Visualization
group_size <- as.numeric(table(cc_obj@idents))

## Aggregated CCC network
save_to <- file.path(plot_dir, paste0(save_token, ".pdf"))
pdf(file = save_to, width = 6, height = 6)
opar <- par(no.readonly = TRUE)
par(mfrow = c(1, 1))
netVisual_circle(cc_obj@net$count, vertex.weight = group_size, weight.scale = T, label.edge= F, title = "Number of interactions")
netVisual_circle(cc_obj@net$weight, vertex.weight = group_size, weight.scale = T, label.edge= F, title = "Interaction weights/strength")
par(opar)
dev.off()

## Per cell type CCC network
save_to <- file.path(plot_dir, paste0(save_token, ".per_celltype.pdf"))
pdf(file = save_to, width = 12, height = 6)
opar <- par(no.readonly = TRUE)
par(mfrow = c(2, 3))
mat <- cc_obj@net$weight
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = group_size, weight.scale = T, edge.weight.max = max(mat))
}
par(opar)
dev.off()

## Circle plot between sender receiver plot, CD8+ T as sender
save_to <- file.path(plot_dir, paste0(save_token, ".circle_plot.pdf"))
pdf(save_to, width = 8, height = 8)
netVisual_chord_gene(cc_obj, sources.use = 3, targets.use = c(1, 2, 4, 5, 6), lab.cex = 0.5, legend.pos.y = 30)
dev.off()

## Bubble plot between sender receiver plot, CD8+ T as sender
save_to <- file.path(plot_dir, paste0(save_token, ".bubble_plot.pdf"))
pdf(save_to, width = 4, height = 3)
netVisual_bubble(cc_obj, sources.use = 3, targets.use = c(1, 2, 4, 5, 6), remove.isolate = FALSE)
dev.off()

## Network centrality scores
save_to <- file.path(plot_dir, paste0(save_token, ".pathway_network.heatmap.pdf"))
pdf(save_to, width = 4, height = 3)
netAnalysis_signalingRole_network(cc_obj)
dev.off()

## Dominate sender receiver plot
save_to <- file.path(plot_dir, paste0(save_token, ".dominate_sender_receiver.dot_plot.pdf"))
pdf(save_to, width = 4, height = 3)
netAnalysis_signalingRole_scatter(cc_obj)
dev.off()

## Signaling contributions
save_to <- file.path(plot_dir, paste0(save_token, ".signaling_contributions.heatmap.pdf"))
pdf(save_to, width = 10, height = 5)
netAnalysis_signalingRole_heatmap(cc_obj, pattern = "outgoing") + netAnalysis_signalingRole_heatmap(cc_obj, pattern = "incoming")
dev.off()

## Estimate communication patterns
save_to <- file.path(plot_dir, paste0(save_token, ".nr_communication_patterns.pdf"))
pdf(save_to, width = 6, height = 3)
selectK(cc_obj, pattern = "outgoing")
dev.off()

# Save the CellChat object
save_to <- file.path(object_dir, paste0(save_token, ".cellchat_obj.rds"))
if (!file.exists(save_to) || overwrite) {
  cat("[I]: Dumping into disk ...\n");
  saveRDS(cc_obj, file = save_to)
}
