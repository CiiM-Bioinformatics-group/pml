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
  library(ComplexHeatmap)
})



rm(list = ls()); gc()

RANDOM_SEED <- 31415926 #
overwrite <- TRUE
proj_dir <- "~/Documents/projects/wp_pml"
plot_dir <- file.path(proj_dir, "outputs/analysis/cell_cell_communication/plots")
object_dir <- file.path(proj_dir, "outputs/analysis/cell_cell_communication/objects")

tar_cell_types <- c("Mono", "CD4 T", "CD8 T", "B", "NK", "DC")

per_omics <- "pml_rnaseq"
per_omics <- "pml_citeseq"

# response <- "Responder"; timepoint <- "BL"
# response <- "Responder"; timepoint <- "6W"
# response <- "Responder"; timepoint <- "3M"
# response <- "Non-responder"; timepoint <- "BL"
# response <- "Non-responder"; timepoint <- "3M"
# save_token <- paste("ccc", per_omics, response, timepoint, sep = ".")


#
## Cell communications by CellChat
#
# CellChat DB of human
cc_db.use <- CellChatDB.human %>% subsetDB(search = "Secreted Signaling")

# Save the CellChat object
sample_table <- tibble::tribble(
  ~Response, ~Timepoint,
  "Responder", "BL",
  "Responder", "3M",
  "Responder", "6W",
  "Non-responder", "BL",
  "Non-responder", "3M",
)
sample_names <- sample_table %>% dplyr::mutate(sample_name = paste(Response, Timepoint, sep = "_")) %>% dplyr::pull(sample_name)
cc_list <- apply(sample_table, 1, function(vec, .per_omics, .tar_cell_types, .overwrite = FALSE) {
  set.seed(RANDOM_SEED)
  response <- vec[1]; timepoint <- vec[2]
  save_token <- paste("ccc", .per_omics, response, timepoint, sep = ".")
  cc_obj_save_to <- file.path(object_dir, paste0(save_token, ".cellchat_obj.rds"))
  if (!file.exists(cc_obj_save_to) || .overwrite) {
    cat("[I]: Creating CellChat object ...", save_token, "\n")
    # Create CellChat object
    pbmc_int <- file.path(proj_dir, "outputs/analysis/integrated/objects", paste0(.per_omics, ".rds")) %>% readRDS()
    tar_cells <- pbmc_int@meta.data %>%
      dplyr::filter(Response == response, Timepoint == timepoint, predicted.celltype.l1 %in% .tar_cell_types) %>%
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

    cat("[I]: Dumping into disk ...\n");
    saveRDS(cc_obj, file = cc_obj_save_to)


    # Visualization
    group_size <- as.numeric(table(cc_obj@idents))

    ## Aggregated CCC network
    save_to <- file.path(plot_dir, paste0(save_token, ".v2.pdf"))
    pdf(file = save_to, width = 6, height = 6)
    netVisual_circle(cc_obj@net$count, vertex.weight = group_size, weight.scale = T, label.edge= F, title = "Number of interactions")
    dev.off()

    save_to <- file.path(plot_dir, paste0(save_token, ".weighted.v2.pdf"))
    pdf(file = save_to, width = 6, height = 6)
    netVisual_circle(cc_obj@net$weight, vertex.weight = group_size, weight.scale = T, label.edge= F, title = "Interaction weights/strength")
    dev.off()

    ## Per cell type CCC network
    save_to <- file.path(plot_dir, paste0(save_token, ".per_celltype.v2.pdf"))
    pdf(file = save_to, width = 12, height = 6)
    opar <- par(no.readonly = TRUE)
    par(mfrow = c(2, 3), mar = c(1, 1, 1, 1))
    mat <- cc_obj@net$weight
    for (i in 1:nrow(mat)) {
      mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
      mat2[i, ] <- mat[i, ]
      netVisual_circle(mat2, vertex.weight = group_size, weight.scale = T, edge.weight.max = max(mat))
    }
    par(opar)
    dev.off()

    ## Circle plot between sender receiver plot, CD8+ T as sender
    save_to <- file.path(plot_dir, paste0(save_token, ".circle_plot.v2.pdf"))
    pdf(save_to, width = 8, height = 8)
    netVisual_chord_gene(cc_obj, sources.use = 3, targets.use = c(1, 2, 4, 5, 6), lab.cex = 0.5, legend.pos.y = 30)
    dev.off()

    ## Bubble plot between sender receiver plot, CD8+ T as sender
    save_to <- file.path(plot_dir, paste0(save_token, ".bubble_plot.v2.pdf"))
    p <- netVisual_bubble(cc_obj, sources.use = 3, targets.use = c(1, 2, 4, 5, 6), remove.isolate = FALSE)
    ggsave(save_to, plot = p, width = 4, height = 3)

    ## Network centrality scores
    save_to <- file.path(plot_dir, paste0(save_token, ".pathway_network.heatmap.v2.pdf"))
    pdf(save_to, width = 4, height = 3)
    netAnalysis_signalingRole_network(cc_obj)
    dev.off()

    ## Dominate sender receiver plot
    save_to <- file.path(plot_dir, paste0(save_token, ".dominate_sender_receiver.dot_plot.pdf"))
    p <- netAnalysis_signalingRole_scatter(cc_obj)
    ggsave(save_to, plot = p, width = 4, height = 3)

    ## Outgoing signaling contributions
    save_to <- file.path(plot_dir, paste0(save_token, ".outgoing_signaling_contributions.heatmap.pdf"))
    pdf(save_to, width = 6, height = 4)
    netAnalysis_signalingRole_heatmap(cc_obj, pattern = "outgoing", width = 6, height = 4) %>% draw()
    dev.off()

    ## Incoming signaling contributions
    save_to <- file.path(plot_dir, paste0(save_token, ".incoming_signaling_contributions.heatmap.pdf"))
    pdf(save_to, width = 6, height = 4)
    netAnalysis_signalingRole_heatmap(cc_obj, pattern = "incoming", width = 6, height = 4) %>% draw()
    dev.off()

    ## Estimate communication patterns
    save_to <- file.path(plot_dir, paste0(save_token, ".nr_communication_patterns.pdf"))
    p <- selectK(cc_obj, pattern = "outgoing")
    ggsave(save_to, plot = p, width = 6, height = 3)
  } else {
    cat("[I]: Loading from disk ...\n");
    cc_obj <- readRDS(cc_obj_save_to)
  }

  cc_obj
}, .per_omics = per_omics, .tar_cell_types = tar_cell_types) %>% purrr::set_names(sample_names)


set.seed(RANDOM_SEED)
# Comparisone between responder and non-responders
# cmp_pairs <- c("NRBL" = "Non-responder_BL", "NR3M" = "Non-responder_3M")
# cmp_pairs <- c("RsBL" = "Responder_BL", "NRBL" = "Non-responder_BL")
# cmp_pairs <- c("Rs6W" = "Responder_6W", "RsBL" = "Responder_BL")
# cmp_pairs <- c("Rs3M" = "Responder_3M", "Rs6W" = "Responder_6W")
# cmp_pairs <- c("Rs3M" = "Responder_3M", "RsBL" = "Responder_BL")

# cmp_pairs <- c("RsBL" = "Responder_BL", "Rs6W" = "Responder_6W", "Rs3M" = "Responder_3M")

save_token <- paste("ccc", per_omics, "comparison", paste0(cmp_pairs, collapse = "_to_"), sep = ".")
cc_obj_save_to <- file.path(object_dir, paste0(save_token, ".cellchat_obj.rds"))
if (file.exists(cc_obj_save_to)) {
  cat("[I]: Loading from disk ...\n");
  cc_obj <- readRDS(cc_obj_save_to)
} else {
  cat("[I]: Dumping to disk ...\n");
  cc_obj <- mergeCellChat(cc_list[cmp_pairs], add.names = names(cmp_pairs))
  saveRDS(cc_obj, file = cc_obj_save_to)
}

# Comparing interaction strength
if (length(cmp_pairs) <= 2) {
  gg1 <- compareInteractions(cc_obj, show.legend = F, group = c(2, 1))
  gg2 <- compareInteractions(cc_obj, show.legend = F, group = c(2, 1), measure = "weight")
  p <- gg1 + gg2
  save_to <- file.path(plot_dir, paste0(save_token, ".communication_stength.pdf"))
  ggsave(save_to, p, width = 6, height = 4)

  # Comparing interaction strength between cell-pairs
  save_to <- file.path(plot_dir, paste0(save_token, ".communication_stength_between_cell_pairs.pdf"))
  pdf(save_to, width = 6, height = 6)
  netVisual_diffInteraction(cc_obj, comparison = c(2, 1), edge.width.max = 15, weight.scale = T)
  dev.off()

  save_to <- file.path(plot_dir, paste0(save_token, ".communication_stength_between_cell_pairs.weighted.pdf"))
  pdf(save_to, width = 6, height = 6)
  netVisual_diffInteraction(cc_obj, comparison = c(2, 1), edge.width.max = 15, weight.scale = T, measure = "weight")
  dev.off()

  # Comparing interaction strength between cell-types
  save_to <- file.path(plot_dir, paste0(save_token, ".communication_stength_between_cell_types.heatmap.pdf"))
  pdf(save_to, width = 6, height = 3)
  netVisual_heatmap(cc_obj, comparison = c(2, 1)) %>% draw()
  dev.off()

  save_to <- file.path(plot_dir, paste0(save_token, ".communication_stength_between_cell_types.heatmap.weighted.pdf"))
  pdf(save_to, width = 6, height = 3)
  netVisual_heatmap(cc_obj, comparison = c(2, 1), measure = "weight") %>% draw()
  dev.off()

  # Identify altered signaling with distinct interaction strength
  save_to <- file.path(plot_dir, paste0(save_token, ".signaling_strength_alteration.pdf"))
  g <- rankNet(cc_obj, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = T, do.stat = TRUE) +
    theme(legend.position = "top")
  ggsave(save_to, g, width = 3, height = 2.5)
}

# Compare LR pairs 
save_to <- file.path(plot_dir, paste0(save_token, ".bubble_plot.v2.pdf"))
p <- netVisual_bubble(cc_obj, sources.use = 3, targets.use = 1:6,  comparison = 1:length(cmp_pairs), angle.x = 45)
ggsave(save_to, plot = p, width = 5.5, height = 4)

# Check the IFNG L+R
cc_obj@meta$datasets <- factor(cc_obj@meta$datasets, levels = names(cmp_pairs)) # set factor level
save_to <- file.path(plot_dir, paste0(save_token, ".IFNG_LR_expression.violin.pdf"))
p <- plotGeneExpression(cc_obj, signaling = "IFN-II", split.by = "datasets", colors.ggplot = T, type = "violin")
ggsave(save_to, plot = p, width = 5, height = 4)

save_to <- file.path(plot_dir, paste0(save_token, ".CCR5_LR_expression.violin.pdf"))
p <- plotGeneExpression(cc_obj, signaling = "CCL", split.by = "datasets", colors.ggplot = T, type = "violin")
ggsave(save_to, plot = p, width = 5, height = 4)
