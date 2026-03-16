#!/usr/bin/env Rscript
# File: deg_analysis.r
# Author: Zhenhua Zhang
# E-mail: zhenhua.zhang217@gmail.com
# Created: Jun 28, 2023
# Updated: Nov 23, 2023

suppressPackageStartupMessages({
  library(gtools)
  library(data.table)
  library(igraph)

  library(tidyverse)
  library(tidygraph)

  library(ggsci)
  library(ggbreak)
  library(ggrepel)
  library(ggVennDiagram)
  library(ggdendro)
  library(ggraph)

  library(circlize)
  library(patchwork)
  library(enrichplot)
  library(RColorBrewer)
  library(ComplexHeatmap)

  library(Seurat)

  library(clusterProfiler)
  library(GOSemSim)
  library(goProfiles)
  library(fgsea)

  library(org.Hs.eg.db)
})

# Utility functions
shrink_to_range <- function(x, t) {
  t_min <- min(t, na.rm = TRUE)
  t_max <- max(t, na.rm = TRUE)
  x_min <- min(x, na.rm = TRUE)
  x_max <- max(x, na.rm = TRUE)
  scaled_x <- t_min + (x - x_min) * (t_max - t_min) / (x_max - x_min)
  return(scaled_x)
}

# colors
rbgc <- brewer.pal(n = 11, name = "RdBu")
# chrom_xy_genes <- c("RPS4Y1", "EIF1AY", "UTY", "EIF1AX", "USP9Y", "ZFY", "XIST", "RPS4X", "JPX", "ZFX", "ZRSR2", "DDX3Y", "DDX3X", "KDM6A", "EIF2S3", "KDM5C", "TXLNG", "PRKX", "CXorf38", "SMC1A", "CA5B", "UBA1", "FUNDC1", "KANTR", "CD99", "SYAP1")
chrom_xy_genes <- fread("~/Documents/projects/wp_pml/inputs/references/sextual_chromosome_genes.txt", header = F) %>% pull(V2)


#
## DEG identification
#
proj_dir <- "~/Documents/projects/wp_pml"
deg_dir <- file.path(proj_dir, "outputs/analysis/deg")
plot_dir <- file.path(deg_dir, "plots")
table_dir <- file.path(deg_dir, "tables")
object_dir <- file.path(proj_dir, "outputs/analysis/integrated/objects")

tar_cell_types <- c("Mono", "CD4 T", "CD8 T", "B", "NK")
celltype_order <- c("Mono", "CD8 T", "CD4 T", "B", "NK", "DC", "other T", "other")

# Colors
col_fun <- colorRamp2(c(0, 1.0), c("gray95", "red"))
hsGO <- godata('org.Hs.eg.db', ont = "BP")

pml_omics_tab <- tibble::tribble(
  ~omics_type, ~fpath,
  "pml_citeseq", "pml_citeseq.rds",
  "pml_rnaseq", "pml_rnaseq.rds",
  "pml_multiome", "pml_multiome.rds",
  # "pml_atacseq", "pml_atacseq.rds", # Not included in the analysis yet.
)

de_tab_save_to <- file.path(proj_dir, "outputs/analysis/deg/tables", paste0("pbmc.integrated.de_gene.csv"))
cp_tab_save_to <- file.path(proj_dir, "outputs/analysis/deg/tables", paste0("pbmc.integrated.cell_proportion.csv"))
if (file.exists(de_tab_save_to, cp_tab_save_to) %>% all()) {
  cat("[I]: Loading from disk... \n")
  de_tab_all <- fread(de_tab_save_to)
  cp_tab_all <- fread(cp_tab_save_to)
} else {
  cat("[I]: Dumping into disk... \n")
  res_list <- apply(pml_omics_tab, 1, function(x, .celltype_order = celltype_order) {
    de_tab <- NULL
    # x <- c('omics_type' = 'pml_citeseq', 'fpath' = 'pml_citeseq.rds')
    omics_type <- x["omics_type"]
    pml_omics <- readRDS(file.path(proj_dir, "outputs/analysis/integrated/objects", x["fpath"]))
    if (omics_type %in% c("pml_citeseq", "pml_multiome")) {
      # Cell proportion
      cell_proportion_tab <- pml_omics@meta.data %>%
        dplyr::group_by(Vireo_assignment, Timepoint, Response, predicted.celltype.l1) %>%
        dplyr::summarise(n = n()) %>%
        dplyr::mutate(prop = n / sum(n)) %>%
        dplyr::arrange(desc(prop)) %>%
        dplyr::mutate(predicted.celltype.l1 = factor(predicted.celltype.l1, levels = .celltype_order)) %>%
        dplyr::mutate(Timepoint = factor(Timepoint, levels = c("BL", "6W", "3M"))) %>%
        dplyr::ungroup() %>%
        dplyr::select(PatientID = Vireo_assignment, Timepoint, Response, predicted.celltype.l1, prop) %>%
        dplyr::filter(!(Timepoint == "3M" & Response == "Non-responder"))

      ## 1. Cell proportion overview
      p <- ggplot(cell_proportion_tab) +
        geom_line(aes(x = Timepoint, y = prop, group = PatientID)) +
        geom_point(aes(x = Timepoint, y = prop, color = Response), alpha = 0.75) +
        scale_fill_npg() +
        theme_classic() +
        theme(axis.text.y = element_text(size = 8), axis.title = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(y = "Proportion", x = "Cell type") +
        facet_grid(~predicted.celltype.l1, space = "free_x")
      save_to <- file.path(proj_dir, "outputs/analysis/cell_proportion/plots", paste0("pbmc.", omics_type, ".integrated.cell_proportion.v2.pdf"))
      ggsave(save_to, plot = p, width = 8.5, height = 3)

      ## 2. Cell proportion by patient
      p <- ggplot(cell_proportion_tab) +
        geom_point(aes(x = PatientID, y = prop, color = predicted.celltype.l1)) +
        scale_color_npg() +
        theme_classic() +
        theme(axis.text.y = element_text(size = 8), axis.title = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(y = "Proportion", x = "Cell type") +
        facet_grid(~Timepoint, scales = "free_x", space = "free_x")
      save_to <- file.path(proj_dir, "outputs/analysis/cell_proportion/plots", paste0("pbmc.", omics_type, ".integrated.cell_proportion_by_patient.pdf"))
      ggsave(save_to, plot = p, width = 8.5, height = 3)

      ## 3. Cell proportion by timepoints
      p <- ggplot(cell_proportion_tab) +
        geom_col(aes(x = Timepoint, y = prop, fill = predicted.celltype.l1), width = 0.985) +
        scale_fill_npg() +
        theme_classic() +
        theme(axis.text.y = element_text(size = 8), axis.title = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(y = "Proportion", x = "Cell type") +
        facet_grid(~PatientID, scales = "free_x", space = "free_x")
      save_to <- file.path(proj_dir, "outputs/analysis/cell_proportion/plots", paste0("pbmc.", omics_type, ".integrated.cell_proportion_by_timepoints.pdf"))
      ggsave(save_to, plot = p, width = 8.5, height = 3)

      # DE gene analysis
      ## 1. DE genes between responder and non-responder per cell type per time point
      Idents(pml_omics) <- paste(pml_omics$predicted.celltype.l1, pml_omics$Response, pml_omics$Timepoint, sep = "_")
      for (pct in tar_cell_types) {
        for (ptp in c("BL", "6W", "3M")) {
          id_1 <- paste0(pct, "_Responder_", ptp)
          id_2 <- paste0(pct, "_Non-responder_", ptp)
          tryCatch({
            de_tab <- FindMarkers(pml_omics, assay = "SCT", ident.1 = id_1, ident.2 = id_2, logfc.threshold = 0, min.pct = 0.05, recorrect_umi = FALSE, verbose = FALSE) %>% 
              dplyr::mutate(gene_symbol = rownames(.), donors = "All", celltype = pct, comparison = paste0("Rs.vs.NRs_", ptp)) %>%
              rbind(de_tab, .)
          }, error = function(e) cat(e$message, "; ", paste(id_1, id_2, "failed, comparison 1\n")))
        }
      }

      ## 2. DE genes between non-responder and responder per cell type per time point. Mering 6W and 3M
      Idents(pml_omics) <- paste(pml_omics$predicted.celltype.l1, pml_omics$Response, ifelse(pml_omics$Timepoint == "BL", "BL", "PT"), sep = "_")
      for (pct in tar_cell_types) {
        id_1 <- paste0(pct, "_Responder_PT")
        id_2 <- paste0(pct, "_Non-responder_PT")
        tryCatch({
            de_tab <- FindMarkers(pml_omics, assay = "SCT", ident.1 = id_1, ident.2 = id_2, logfc.threshold = 0, min.pct = 0.05, recorrect_umi = FALSE, verbose = FALSE) %>%
              dplyr::mutate(gene_symbol = rownames(.), donors = "All", celltype = pct, comparison = "Rs.vs.NRs_PT") %>%
              rbind(de_tab)
        }, error = function(e) cat(e$message, "; ", paste(id_1, id_2, "failed, comparison 2\n")))
      }

      ## 3. DE genes between three time points in responders per cell type.
      Idents(pml_omics) <- paste(pml_omics$predicted.celltype.l1, pml_omics$Timepoint, sep = "_")
      tar_cells <- pml_omics@meta.data %>% dplyr::filter(Response == "Responder") %>% rownames()
      for (pct in tar_cell_types) {
        id_1 <- paste0(pct, "_BL")
        id_2 <- paste0(pct, "_6W")
        id_3 <- paste0(pct, "_3M")
        tryCatch({
          de_tab <- FindMarkers(pml_omics[, tar_cells], assay = "SCT", ident.1 = id_3, ident.2 = id_1, logfc.threshold = 0, min.pct = 0.05, recorrect_umi = FALSE, verbose = FALSE) %>%
            dplyr::mutate(gene_symbol = rownames(.), donors = "Responder", celltype = pct, comparison = "3M.vs.BL") %>%
            rbind(de_tab)
        }, error = function(e) cat(e$message, "; ", paste(id_2, id_1, "failed, comparison 3\n")))
        tryCatch({
          de_tab <- FindMarkers(pml_omics[, tar_cells], assay = "SCT", ident.1 = id_2, ident.2 = id_1, logfc.threshold = 0, min.pct = 0.05, recorrect_umi = FALSE, verbose = FALSE) %>%
            dplyr::mutate(gene_symbol = rownames(.), donors = "Responder", celltype = pct, comparison = "6W.vs.BL") %>%
            rbind(de_tab, .)
        }, error = function(e) cat(e$message, "; ", paste(id_2, id_1, "failed, comparison 3\n")))
        tryCatch({
          de_tab <- FindMarkers(pml_omics[, tar_cells], assay = "SCT", ident.1 = id_3, ident.2 = id_2, logfc.threshold = 0, min.pct = 0.05, recorrect_umi = FALSE, verbose = FALSE) %>%
            dplyr::mutate(gene_symbol = rownames(.), donors = "Responder", celltype = pct, comparison = "3M.vs.6W") %>%
            rbind(de_tab, .)
        }, error = function(e) cat(e$message, "; ", paste(id_2, id_1, "failed, comparison 3\n")))
      }

      ## 4. DE genes between three time points in non-responders per cell type.
      Idents(pml_omics) <- paste(pml_omics$predicted.celltype.l1, pml_omics$Timepoint, sep = "_")
      tar_cells <- pml_omics@meta.data %>% dplyr::filter(Response == "Non-responder") %>% rownames()
      for (pct in tar_cell_types) {
        id_1 <- paste0(pct, "_BL")
        id_2 <- paste0(pct, "_3M")
        tryCatch({
          de_tab <- FindMarkers(pml_omics[, tar_cells], assay = "SCT", ident.1 = id_2, ident.2 = id_1, logfc.threshold = 0, min.pct = 0.05, recorrect_umi = FALSE, verbose = FALSE) %>%
            dplyr::mutate(gene_symbol = rownames(.), donors = "Non-responder", celltype = pct, comparison = "3M.vs.BL") %>%
            rbind(de_tab, .)
        }, error = function(e) cat(e$message, "; ", paste(id_2, id_1, "failed, comparison 4\n")) )
      }
    } else if (omics_type == "pml_rnaseq") {
      cell_proportion_tab <- pml_omics@meta.data %>%
        dplyr::group_by(Vireo_assignment, Timepoint, Response, predicted.celltype.l1) %>%
        dplyr::summarise(n = n()) %>%
        dplyr::mutate(prop = n / sum(n)) %>%
        dplyr::arrange(desc(prop)) %>%
        dplyr::mutate(predicted.celltype.l1 = factor(predicted.celltype.l1, levels = .celltype_order)) %>%
        dplyr::mutate(Timepoint = factor(Timepoint, levels = c("BL", "6W", "3M"))) %>%
        dplyr::ungroup() %>%
        dplyr::select(PatientID = Vireo_assignment, Timepoint, Response, predicted.celltype.l1, prop)

      # 1. Cell proportion overview
      p <- ggplot(cell_proportion_tab) +
        geom_line(aes(x = Timepoint, y = prop, group = PatientID)) +
        geom_point(aes(x = Timepoint, y = prop, color = Response), alpha = 0.75) +
        scale_fill_npg() +
        theme_classic() +
        theme(axis.text.y = element_text(size = 8), axis.title = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(y = "Proportion", x = "Cell type") +
        facet_grid(~predicted.celltype.l1, space = "free_x")
      save_to <- file.path(proj_dir, "outputs/analysis/cell_proportion/plots", paste0("pbmc.", omics_type, ".integrated.cell_proportion.pdf"))
      ggsave(save_to, plot = p, width = 8.5, height = 3)

      # 2. Cell proportion by patient
      p <- ggplot(cell_proportion_tab) +
        geom_point(aes(x = PatientID, y = prop, color = predicted.celltype.l1)) +
        scale_color_npg() +
        theme_classic() +
        theme(axis.text.y = element_text(size = 8), axis.title = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(y = "Proportion", x = "Cell type") +
        facet_grid(~Timepoint, scales = "free_x", space = "free_x")
      save_to <- file.path(proj_dir, "outputs/analysis/cell_proportion/plots", paste0("pbmc.", omics_type, ".integrated.cell_proportion_by_patient.pdf"))
      ggsave(save_to, plot = p, width = 8.5, height = 3)

      # 3. Cell proportion by timepoints
      p <- ggplot(cell_proportion_tab) +
        geom_col(aes(x = Timepoint, y = prop, fill = predicted.celltype.l1), width = 0.985) +
        scale_fill_npg() +
        theme_classic() +
        theme(axis.text.y = element_text(size = 8), axis.title = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(y = "Proportion", x = "Cell type") +
        facet_grid(~PatientID, scales = "free_x", space = "free_x")
      save_to <- file.path(proj_dir, "outputs/analysis/cell_proportion/plots", paste0("pbmc.", omics_type, ".integrated.cell_proportion_by_timepoints.pdf"))
      ggsave(save_to, plot = p, width = 7, height = 3)

      # DE gene analysis
      ## 1. DE genes between responder and non-responder per cell type per time point
      Idents(pml_omics) <- paste(pml_omics$predicted.celltype.l1, pml_omics$Response, pml_omics$Timepoint, sep = "_")
      for (pct in tar_cell_types) {
        for (ptp in c("BL")) {
          id_1 <- paste0(pct, "_Responder_", ptp)
          id_2 <- paste0(pct, "_Non-responder_", ptp)
          tryCatch({
            de_tab <- FindMarkers(pml_omics, ident.1 = id_1, ident.2 = id_2, assay = "SCT", logfc.threshold = 0, min.pct = 0.05, recorrect_umi = FALSE, verbose = FALSE) %>% 
              dplyr::mutate(gene_symbol = rownames(.), donors = "All", celltype = pct, comparison = paste0("Rs.vs.NRs_", ptp)) %>%
              rbind(de_tab, .)
          }, error = function(e) cat(e$message, "; ", paste(id_1, id_2, "failed, comparison 1\n")))
        }
      }
    }

    cell_proportion_tab$data_source <- omics_type
    save_to <- file.path(proj_dir, "outputs/analysis/cell_proportion/tables", paste0("pbmc.", omics_type, "integrated.cell_propotion.csv"))
    fwrite(cell_proportion_tab, save_to)

    de_tab <- dplyr::mutate(de_tab, direction = dplyr::if_else(avg_log2FC > 0, "up", "dw"), data_source = omics_type)
    save_to <- file.path(proj_dir, "outputs/analysis/deg/tables", paste0("pbmc.", omics_type, ".integrated.de_gene.csv"))
    fwrite(de_tab, save_to)

    list(detab = de_tab, cptab = cell_proportion_tab)
  }, .celltype_order = celltype_order)


  ## Collection all DE results
  de_tab_all <- res_list %>% lapply(function(x) x$detab) %>% Reduce(rbind, .)
  fwrite(de_tab_all, de_tab_save_to)
  cp_tab_all <- res_list %>% lapply(function(x) x$cptab) %>% Reduce(rbind, .)
  fwrite(cp_tab_all, cp_tab_save_to)
}


# 1. Plot number of DE genes between responders and non-responders at each time point per cell type.
per_omics <- "pml_rnaseq"
per_omics <- "pml_citeseq"
de_fmt_tab <- as.data.table(de_tab_all) %>%
  dplyr::filter(p_val_adj <= 0.05, abs(avg_log2FC) >= 0.1, pct.2 >= 0.05, pct.1 >= 0.05, donors == "All", data_source == per_omics, !(gene_symbol %in% chrom_xy_genes)) %>% 
  tidyr::separate(comparison, into = c("comparison", "timepoint"), sep = "_") %>%
  dplyr::mutate(timepoint = factor(timepoint, levels = c("BL", "3M", "PT")))


# Bar plot show number of DE genes
de_count_tab <- de_fmt_tab %>%
  dplyr::filter(data_source == per_omics, timepoint != "PT") %>%
  dplyr::group_by(data_source, celltype, timepoint, comparison, direction) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(celltype = forcats::fct_reorder(celltype, -n)) %>%
  dplyr::mutate(direction = factor(direction, levels = c("up", "dw")))

p <- de_count_tab %>%
  dplyr::filter(data_source == per_omics) %>%
  ggplot() +
  geom_col(aes(x = celltype, y = n, fill = direction), width = 0.8, position = position_dodge2(width = 0.75, preserve = "single")) +
  geom_text(aes(x = celltype, y = n, group = direction, label = n, vjust = -0.1), size = 4, position = position_dodge2(width = 0.75, preserve = "single")) +
  scale_fill_npg() +
  labs(x = NULL, y = NULL, fill = "Responders vs non-responders") +
  theme_classic() +
  theme(legend.position = "top", axis.text.y = element_text(size = 8), axis.text.x = element_text(angle = 45, hjust = 1), strip.background = element_blank()) +
  facet_wrap(~timepoint, nrow = 1)
save_to <- file.path(proj_dir, "outputs/analysis/deg/plots", paste0("de_gene.", per_omics, ".Rs_vs_NR.per_timepoint.pdf"))
ggsave(save_to, plot = p, width = 8, height = 4)

# Venn diagram to show DE genes shared
for (pct in tar_cell_types) {
  venn_plot_list <- de_fmt_tab %>%
    dplyr::filter(timepoint %in% c("BL", "3M"), celltype == pct) %>%
    dplyr::select(gene_symbol, direction, timepoint) %>%
    (function(tab) {
      list(
        `BL (up)` = tab %>% dplyr::filter(direction == "up", timepoint == "BL") %>% dplyr::pull(gene_symbol),
        `3M (up)` = tab %>% dplyr::filter(direction == "up", timepoint == "3M") %>% dplyr::pull(gene_symbol),
        `BL (dw)` = tab %>% dplyr::filter(direction == "dw", timepoint == "BL") %>% dplyr::pull(gene_symbol),
        `3M (dw)` = tab %>% dplyr::filter(direction == "dw", timepoint == "3M") %>% dplyr::pull(gene_symbol)
      )})

  pct <- stringr::str_remove_all(pct, " ")
  save_to <- file.path(proj_dir, "outputs/analysis/deg/plots", paste("de_gene", per_omics, "per_timepoint.venn_diagram.v2.pdf", sep = "."))
  p <- ggVennDiagram(venn_plot_list) + scale_fill_gradient(low = "white", high = "darkblue") + theme(legend.position = "none")
  ggsave(save_to, plot = p, width = 6, height = 6)
}

# Overlaps of identified DEGs
comb_mat <- de_fmt_tab %>%
    dplyr::filter(timepoint == "BL") %>%
    dplyr::select(gene_symbol, direction, celltype) %>%
    dplyr::group_by(celltype, direction) %>%
    dplyr::summarise(gene_sets =  list(gene_symbol)) %>%
    dplyr::mutate(gene_set_name = paste0(stringr::str_remove_all(celltype, " "), " (", direction, ")")) %>%
    dplyr::pull(gene_sets, gene_set_name) %>%
    list_to_matrix() %>%
    make_comb_mat()

gene_upset_save_to <- file.path(plot_dir, paste("de_gene", per_omics, "all_celltype.Rs_vs_NR.BL.upset_plot.pdf", sep = ".")) %>% stringr::str_remove_all(" ")
pdf(gene_upset_save_to, width = 6.5, height = 3.25)
UpSet(comb_mat,top_annotation = upset_top_annotation(comb_mat, add_numbers = TRUE), left_annotation = upset_left_annotation(comb_mat, add_numbers = TRUE), bg_col = "#F0F0FF", bg_pt_col = "#CCCCFF") %>% draw()
dev.off()


# Volcano plot of DEGs of each cell type
for (pct in tar_cell_types) {
  de_tab_per_ct <- de_tab_all %>%
    dplyr::filter(data_source == "pml_citeseq", celltype == pct, donors == "All") %>%
    tidyr::separate(comparison, into = c("comparison", "timepoint"), sep = "_") %>%
    dplyr::filter(timepoint %in% c("BL"), ! gene_symbol %in% chrom_xy_genes) %>%
    dplyr::select(p_val_adj, avg_log2FC, gene_symbol) %>%
    dplyr::mutate(Signif. = dplyr::case_when(p_val_adj <= 0.05 & avg_log2FC <= -0.25 ~ "Down", p_val_adj <= 0.05 & avg_log2FC >= 0.25 ~ "Up", TRUE ~ "Other"))

  volp <- ggplot() +
    geom_point(aes(x = avg_log2FC, y = -log10(p_val_adj), color = Signif.), de_tab_per_ct) +
    geom_text_repel(aes(x = avg_log2FC, y = -log10(p_val_adj), label = gene_symbol), dplyr::slice_min(de_tab_per_ct, p_val_adj, n = 20), max.overlaps = Inf, min.segment.length = 0) +
    geom_vline(xintercept = c(0.25, -0.25), linetype = "dotted") +
    geom_hline(yintercept = -log10(0.05), linetype = "dotted") +
    labs(x = "Log2(fold-change)", y = "-log10(Adj. P value)") +
    scale_color_manual(values = c(Up = "red", Down = "blue", Other = "grey")) +
    theme_classic() +
    theme(legend.position = "top")

  volp_save_to <- file.path(proj_dir, "outputs/analysis/deg/plots", paste("de_gene", per_omics, stringr::str_remove_all(pct, " "), "Rs_vs_NR.BL.volcano_plot.pdf", sep = "."))
  ggsave(volp_save_to, plot = volp, width = 4, height = 5)
}


# DE genes replicated by RNA-seq
plot_tab <- de_tab_all %>%
  dplyr::filter(
    donors == "All", comparison == "Rs.vs.NRs_BL",
    data_source %in% c("pml_citeseq", "pml_rnaseq"),
    pct.1 >= 0.05, pct.2 >= 0.05,
    abs(avg_log2FC) > 0.01
  ) %>%
  dplyr::select(data_source, celltype, gene_symbol, direction, avg_log2FC, p_val, p_val_adj) %>%
  tidyr::pivot_wider(id_cols = c(celltype, gene_symbol), names_from = c(data_source), values_from = c(direction, avg_log2FC, p_val, p_val_adj)) %>%
  dplyr::mutate(Replication = dplyr::if_else(direction_pml_citeseq == direction_pml_rnaseq, "Replicated", "Not replicated")) %>%
  dplyr::filter(!is.na(Replication))

p <- ggplot(plot_tab) +
  geom_point(aes(x = avg_log2FC_pml_citeseq, y = avg_log2FC_pml_rnaseq, color = Replication)) +
  facet_wrap(~celltype, scales = "free") +
  theme_classic()
save_to <- file.path(plot_dir, "de_gene.replicated_at_BL.points.pdf")
ggsave(save_to, plot = p, width = 10, height = 10)


dis_deg_tab <- de_tab_all %>%
  dplyr::filter(donors == "All", comparison == "Rs.vs.NRs_BL", data_source == "pml_citeseq", p_val_adj < 0.05, abs(avg_log2FC) > 0.1, pct.1 > 0.05, pct.2 > 0.05) %>%
  dplyr::select(celltype, gene_symbol, direction, p_val, avg_log2FC)

val_deg_tab <- de_tab_all %>%
  dplyr::filter(donors == "All", comparison == "Rs.vs.NRs_BL", data_source == "pml_rnaseq") %>%
  dplyr::filter(p_val < 0.05) %>%
  dplyr::select(celltype, gene_symbol, direction, p_val, avg_log2FC)

rep_count_tab <- dis_deg_tab %>%
  dplyr::left_join(val_deg_tab, by = c("celltype", "gene_symbol", "direction"), suffix = c(".citeseq", ".rnaseq")) %>%
  dplyr::mutate(is_replicated = (!is.na(p_val.rnaseq)) & (avg_log2FC.citeseq * avg_log2FC.rnaseq > 0)) %>%
  dplyr::group_by(celltype, direction) %>%
  dplyr::summarise(n = n(), n_rep = sum(is_replicated)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(celltype = forcats::fct_reorder(celltype, -n)) %>%
  dplyr::mutate(direction = factor(direction, levels = c("up", "dw")), timepoint = "BL")

p <- ggplot(rep_count_tab) +
  geom_col(aes(x = celltype, y = n, group = direction), fill = "gray", width = 0.8, position = position_dodge2(width = 0.75, preserve = "single")) +
  geom_col(aes(x = celltype, y = n_rep, fill = direction), width = 0.8, position = position_dodge2(width = 0.75, preserve = "single")) +
  geom_text(aes(x = celltype, y = n_rep, label = n_rep, group = direction), position = position_dodge2(width = 0.75, preserve = "single"), vjust = -0.2) +
  geom_text(aes(x = celltype, y = n, label = n, group = direction), position = position_dodge2(width = 0.75, preserve = "single"), vjust = -0.2) +
  labs(x = NULL, y = NULL, fill = "Responder vs non-responder") +
  scale_fill_npg() +
  theme_classic() +
  facet_wrap(~timepoint, nrow = 1, strip.position = "right") +
  theme(legend.position = "top", axis.text.y = element_text(size = 8), axis.text.x = element_text(angle = 45, hjust = 1), strip.background = element_blank())
save_to <- file.path(plot_dir, "de_gene.replicated_at_BL.pdf")
ggsave(save_to, plot = p, width = 4, height = 4.5)


# DEG concordance dynamics compared to Rs vs NR at BL
for (pct in c("CD8 T", "CD4 T", "NK", "Mono", "B")) {
  deg_cordance_tab <- de_tab_all %>%
    dplyr::filter(celltype == pct, data_source == "pml_citeseq", comparison %in% c("3M.vs.6W", "6W.vs.BL", "Rs.vs.NRs_BL"), donors %in% c("All", "Responder"), pct.1 > 0.1, pct.2 > 0.1) %>%
    dplyr::select(gene_symbol, p_val_adj, avg_log2FC, comparison)
  deg_orders <- deg_cordance_tab %>% dplyr::filter(comparison == "Rs.vs.NRs_BL", !gene_symbol %in% chrom_xy_genes) %>% dplyr::arrange(avg_log2FC) %>% dplyr::pull(gene_symbol)
  deg_cordance_plot_tab <- deg_cordance_tab %>% dplyr::filter(gene_symbol %in% deg_orders) %>% dplyr::mutate(gene_symbol = factor(gene_symbol, levels = deg_orders), comparison = factor(comparison, c("Rs.vs.NRs_BL", "6W.vs.BL", "3M.vs.6W")))

  hl_tab <- dplyr::filter(deg_cordance_plot_tab, p_val_adj < 0.05, abs(avg_log2FC) > 0.1)
  p_scatter <- deg_cordance_plot_tab %>% ggplot() +
    geom_hline(yintercept = 0, linewidth = 0.25) +
    geom_point(aes(x = gene_symbol, y = avg_log2FC), color = "grey", size = 0.25, alpha = 0.5) +
    geom_point(aes(x = gene_symbol, y = avg_log2FC), dplyr::filter(hl_tab, comparison == "Rs.vs.NRs_BL"), color = "darkred", size = 0.5, alpha = 0.75) +
    geom_point(aes(x = gene_symbol, y = avg_log2FC), dplyr::filter(hl_tab, comparison == "6W.vs.BL"), color = "darkred", size = 0.5, alpha = 0.75) +
    geom_point(aes(x = gene_symbol, y = avg_log2FC), dplyr::filter(hl_tab, comparison == "3M.vs.6W"), color = "darkred", size = 0.5, alpha = 0.75) +
    scale_x_discrete(expand = c(unit(0.025, "npc"), unit(0.025, "npc"))) +
    scale_color_manual(values = c("Y" = "darkred", "N" = "grey")) +
    facet_grid(comparison ~ .) +
    labs(x = "Genes order by log2FC", y = "Log2FC") +
    theme_classic() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  p_scater_saveto <- file.path(plot_dir, paste0("de_gene.pml_citeseq.", stringr::str_remove(pct, " "), ".concordance_across_comparisons.pdf"))
  ggsave(p_scater_saveto, plot = p_scatter, width = 7, height = 4)
}


# Check the direction of DE genes
# Expression dynamic patterns
pct <- "CD8 T"

pbmc_int <- readRDS(file.path(proj_dir, "outputs/analysis/integrated/objects/pml_citeseq.rds"))
tar_features <- deg_cordance_plot_tab %>% dplyr::filter(p_val_adj < 0.05, comparison == "Rs.vs.NRs_BL", abs(avg_log2FC) > 0.2) %>% dplyr::pull(gene_symbol)
norm_avg_expr_tab <- AverageExpression(pbmc_int, assay = "SCT", features = tar_features, group.by = c("Response", "Timepoint"))[["SCT"]] %>%
  as.data.frame() %>%
  tibble::rownames_to_column("gene_symbol") %>%
  dplyr::mutate(expr_pattern = apply(., 1, function(x) paste0(ifelse(x["Non-responder_BL"] < x["Responder_BL"], "+", "-"), ifelse(x["Responder_BL"] < x["Responder_6W"], "+", "-"), ifelse(x["Responder_6W"] < x["Responder_3M"], "+", "-")))) %>%
  tidyr::pivot_longer(-c(gene_symbol, expr_pattern), names_to = "SampleGroup", values_to = "AverageExpression") %>%
  dplyr::filter(!SampleGroup %in% c("Non-responder_3M"), !gene_symbol %in% c("MALAT1", "FTH1", "FTL")) %>%
  dplyr::mutate(AverageExpression = log10(AverageExpression + 1) %>% scale() %>% as.numeric()) %>%
  dplyr::mutate(SampleGroup = factor(SampleGroup, levels = c("Non-responder_BL", "Responder_BL", "Responder_6W", "Responder_3M"))) %>%
  dplyr::mutate(expr_pattern = factor(expr_pattern, levels = c("+++", "++-", "---", "-++", "+--", "--+", "-+-", "+-+")))

hc_method <- "ward.D2"; dst_method <- "euclidean"; ct_k <- 15
gene_symbol_hc <- norm_avg_expr_tab %>%
  tidyr::pivot_wider(names_from = "SampleGroup", values_from = "AverageExpression") %>%
  tibble::column_to_rownames("gene_symbol") %>%
  as.matrix() %>%
  dist(method = dst_method) %>%
  hclust(method = hc_method)
gene_order <- gene_symbol_hc$labels[gene_symbol_hc$order]

norm_avg_expr_plot_tab <- norm_avg_expr_tab %>% dplyr::mutate(gene_symbol = factor(gene_symbol, gene_order))
p <- ggplot(norm_avg_expr_plot_tab) +
  geom_tile(aes(x = SampleGroup, y = gene_symbol, fill = AverageExpression)) +
  labs(x = NULL, y = "Gene symbol", fill = "Avg. expression") +
  scale_fill_gradient2(low = rbgc[11], mid = rbgc[6], high = rbgc[1]) +
  facet_grid(expr_pattern ~ ., scales = "free_y", space = "free_y",) +
  theme_classic() +
  theme(legend.position = "top", panel.spacing = unit(0.002, "npc"), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.text.x = element_text(angle = -30, hjust = 0, vjust = 1))
plot_save_to <- file.path(plot_dir, paste0("de_gene.pml_citeseq.", stringr::str_remove(pct, " "), ".dynamics_patterns.pdf"))
ggsave(plot_save_to, plot = p, width = 4.5, height = 9.75)

# tar_cells <- pbmc_int@meta.data %>% as.data.frame() %>% dplyr::filter(Response != "Non-responder" | (Response == "Non-responder" & Timepoint == "BL"), predicted.celltype.l1 == pct) %>% rownames()
# Idents(pbmc_int) <- paste(pbmc_int$Response, pbmc_int$Timepoint, sep = "_")
# p <- DotPlot(pbmc_int[, tar_cells], features = tar_features, assay = "SCT") +
#   scale_y_discrete(limits = c("Non-responder_BL", "Responder_BL", "Responder_6W", "Responder_3M")) +
#   scale_color_gradient2(low = rbgc[11], mid = rbgc[6], high = rbgc[1]) +
#   coord_flip() +
#   RotatedAxis()
# example_feature_plot_save_to <- file.path(plot_dir, paste0("de_gene.pml_citeseq.", stringr::str_remove(pct, " "), ".concordance_across_comparisons.example_features.same_direction.pdf"))
# ggsave(example_feature_plot_save_to, plot = p, width = 5, height = 12)


# GSEA analysis
ref_de_genes <- de_tab_all %>%
  dplyr::filter(data_source == "pml_citeseq", donors == "All", p_val_adj < 0.1, abs(avg_log2FC) > 0.1, ! gene_symbol %in% chrom_xy_genes) %>%
  dplyr::group_by(celltype) %>%
  dplyr::summarise(gene_set = list(gene_symbol)) %>%
  dplyr::mutate(gene_set_name = stringr::str_remove_all(celltype, " ")) %>%
  dplyr::pull(gene_set, gene_set_name)

for (pct in c("Mono", "CD4 T", "CD8 T", "NK", "B")) {
  for (per_cmp in c("3M.vs.BL", "6W.vs.BL", "3M.vs.6W")) {
    dym_de_genes <- de_tab_all %>%
      dplyr::filter(data_source == "pml_citeseq", donors == "Responder", comparison == per_cmp, p_val_adj < 0.05, abs(avg_log2FC) > 0.1, pct.1 > 0.05, pct.2 > 0.05, celltype == pct, ! gene_symbol %in% chrom_xy_genes) %>%
      dplyr::pull(avg_log2FC, gene_symbol) %>%
      sort(decreasing = TRUE)

    per_gsea_res <- fgsea(ref_de_genes, dym_de_genes, minSize = 10, maxSize = 1000, eps = 0.0)
    gsea_table_save_to <- file.path(plot_dir, "GSEA", paste0("deg.", stringr::str_remove_all(pct, " "), ".", per_cmp, "_dym_deg_in_BL.pdf"))
    pdf(gsea_table_save_to, width = 7, height = 3.5)
    plotGseaTable(ref_de_genes, dym_de_genes, per_gsea_res, colwidths = c(2, 3, 0.8, 1.2, 1.2), gseaParam = 0.5)
    dev.off()
  }
}


# Gene expression alteration pattern
exp_alt_pat_save_to <- file.path(plot_dir, paste0("de_gene.", per_omics, ".all_celltype.Rs_vs_NR.BL.expression_alteration_pattern.v3.pdf"))
if (!file.exists(exp_alt_pat_save_to)) {
  cat("[I]: File exists ... \n")
} else {
  hc_method <- "ward.D2"; dst_method <- "euclidean"; ct_k <- 10
  selected_features <- de_tab_all %>%
    dplyr::filter(!(gene_symbol %in% chrom_xy_genes), p_val_adj <= 0.05, abs(avg_log2FC) >= 0.3, pct.2 >= 0.05, pct.1 >= 0.05, donors == "All", data_source == per_omics) %>%
    dplyr::filter(comparison %in% "Rs.vs.NRs_BL") %>%
    dplyr::pull(gene_symbol) %>%
    unique()

  exp_alt_pat_tab <- de_tab_all %>%
    dplyr::filter(gene_symbol %in% selected_features, donors == "All", data_source == per_omics, comparison == "Rs.vs.NRs_BL") %>%
    dplyr::select(gene_symbol, celltype, p_val_adj, avg_log2FC, direction)

  plot_tab <- exp_alt_pat_tab %>%
    dplyr::select(gene_symbol, celltype, avg_log2FC) %>%
    tidyr::pivot_wider(names_from = "celltype", values_from = "avg_log2FC") %>%
    dplyr::mutate(dplyr::across(-gene_symbol, ~dplyr::if_else(is.na(.x), 0, .x))) %>%
    tibble::column_to_rownames("gene_symbol") %>%
    as.matrix()

  set.seed(31415)
  idx_five_mf <- c("CXCL2", "CXCL3", "CXCL8", "IFITM3", "IFI44L", "MARCH1")
  idx_seven_mf <- c("CD52", "CD8A", "CD8B", "CD226",  "ITGA4", "GZMB", "AKT3")
  idx_eight_mf <- c("IFNG", "IFITM2", "LAIR2", "TGFBR3", "ITK", "TNF", "PRMT2", "HOPX", "ZEB2")
  idx_ten_mf <- c("GNLY", "GZMH", "ITGB1", "CD2", "CCL4")
  mark_features <- c(idx_five_mf, idx_seven_mf, idx_eight_mf, idx_ten_mf)
  mark_colors <- c(
    rep("darkgreen", length(idx_five_mf)), rep("darkred", length(idx_seven_mf)),
    rep("darkblue", length(idx_eight_mf)), rep("orange", length(idx_ten_mf)))
  mark_idx <- match(mark_features, rownames(plot_tab))

  col_fun <- colorRamp2(c(-1, 0, 1), c("#3C5488B2", "white", "#E64B35B2"))

  panel_colors <- c("gray", "gray", "gray", "gray", "darkgreen", "darkblue", "darkred", "gray", "gray", "orange")
  left_annotation <- rowAnnotation(
    f = anno_block(gp = gpar(fill = panel_colors), labels_gp = gpar(col = "white", fontsize = 10))
  )
  right_annotation <- rowAnnotation(
    b = anno_mark(at = mark_idx, labels = mark_features, labels_gp = gpar(fontsize = 10, fontface = "italic", col = mark_colors))
  )

  pdf(exp_alt_pat_save_to, width = 5, height = 11)
  heat_map <- Heatmap(
    plot_tab, col = col_fun, name = "log2FC",
    cluster_columns = FALSE, cluster_rows = TRUE, show_row_names = FALSE, row_km = 10,
    left_annotation = left_annotation, right_annotation = right_annotation
  )
  heat_map <- heat_map %>% draw()
  rcl.list <- row_order(heat_map)
  dev.off()
}


# Functional enrichment, gene ontology enrichment
go_tab_save_to <- file.path(table_dir, paste0("de_gene.all_celltype.go_enrichment.csv"))
if (file.exists(go_tab_save_to)) {
  cat("[W]: Loading from disk ...", go_tab_save_to, "\n")
  go_tab <- fread(go_tab_save_to)
} else {
  de_gene_list <- de_tab_all %>% dplyr::filter(p_val_adj <= 0.05, abs(avg_log2FC) >= 0.1, pct.2 >= 0.05, pct.1 >= 0.05, !(gene_symbol %in% chrom_xy_genes)) %>% 
    dplyr::select(data_source, celltype, comparison, direction, gene_symbol) %>%
    dplyr::group_by(data_source, celltype, comparison, direction) %>%
    dplyr::summarise(gene_list = list(gene_symbol)) %>%
    dplyr::mutate(gene_list_name = paste(data_source, celltype, direction, comparison, sep = "_") %>% stringr::str_remove_all("pml_| ")) %>%
    dplyr::pull(gene_list, gene_list_name)

  go_tab <- lapply(names(de_gene_list), function(x, .gsl) {
    gene_vec <- .gsl[[x]]
    tryCatch({
      ego <- enrichGO(gene_vec, org.Hs.eg.db, "SYMBOL", "ALL")
      dplyr::mutate(ego@result, gene_set = x, n_genes = length(gene_vec)) %>%
        tidyr::separate(gene_set, into = c("data_source", "celltype", "direction", "comparison"), sep = "_", extra = "merge") %>%
        dplyr::mutate(data_source = paste0("pml_", data_source))
    }, error = function(e) { cat(e$message, paste(x, "failed\n")); NULL })
  }, .gsl = de_gene_list) %>%
    Reduce(rbind, .)

  fwrite(go_tab, go_tab_save_to)
}

# Gene ontology enrichment and similarity
# 'arg' should be one of “Wang”, “Resnik”, “Rel”, “Jiang”, “Lin”, “TCSS”
go_sim_mat <- go_tab %>%
  dplyr::filter(data_source == per_omics, comparison == "Rs.vs.NRs_BL", ONTOLOGY == "BP", p.adjust <= 0.05, Count >= 3) %>%
  dplyr::pull(ID) %>%
  as.character() %>%
  unique() %>%
  mgoSim(GO1 = ., GO2 = ., semData = hsGO, measure = "Rel", combine = NULL) %>%
  (function(mat) { mat[is.na(mat)] <- 0; mat })

col_fun <- colorRamp2(c(0, 0.75), c("#90D5FFCC", "darkblue"))
go_sim_p_saveto <- file.path(plot_dir, paste0("de_gene.", per_omics, ".all_celltype.BP.go_sim_tree.v3.pdf"))
pdf(go_sim_p_saveto, width = 8, height = 7)
heat_map <- Heatmap(
  go_sim_mat, col = col_fun, name = "GO similarity",
  cluster_rows = TRUE, show_row_names = FALSE,# row_km = 3,
  cluster_columns = TRUE, show_column_names = FALSE#, column_km = 3
  #left_annotation = left_annotation, right_annotation = right_annotation
) %>% draw()
dev.off()


dist_method <- "manhattan"; cls_method <- "ward.D2"
go_dist <- dist(go_sim_mat, dist_method)
go_cluster <- hclust(go_dist, cls_method)
go_term_sim_order <- go_cluster$labels[go_cluster$order]
go_sim_mat <- go_sim_mat %>% as.data.frame() %>%
  dplyr::mutate(ID1 = rownames(.)) %>%
  tidyr::pivot_longer(cols = -ID1, names_to = "ID2", values_to = "similarity") %>%
  dplyr::mutate(similarity = dplyr::if_else(ID1 == ID2, 1.0, similarity)) %>%
  dplyr::mutate(ID1 = factor(ID1, levels = go_term_sim_order)) %>%
  dplyr::mutate(ID2 = factor(ID2, levels = go_term_sim_order)) %>%
  dplyr::mutate(ID_order = as.integer(ID1))

go_term_order_tab <- go_sim_mat %>% dplyr::select(GO_ID = ID1, GO_ID_index = ID_order) %>% dplyr::distinct()

ct_k <- 5
go_cluster_tree <- cutree(go_cluster, k = ct_k)
go_cluster_hc <- as.dendrogram(go_cluster)
go_cluster_dd <- dendro_data(go_cluster_hc, type = "rectangle")

go_example_hc <- data.frame(cluster = factor(go_cluster_tree)) %>% tibble::rownames_to_column("label")
go_cluster_dd[["labels"]] <- merge(go_cluster_dd[["labels"]], go_example_hc, by="label")
go_cluster_tree_rect <- go_cluster_dd[["labels"]] %>% dplyr::group_by(cluster) %>% dplyr::summarise(Y.1 = min(x) - 0.5, Y.2 = max(x) + 0.5)
go_cluster_tree_xmax <- mean(go_cluster$height[length(go_cluster$height)-((ct_k-2):(ct_k-1))]) %>% sqrt()

go_pat_tab <- selected_go_tab %>%
  dplyr::select(ID, celltype, direction, Count) %>%
  dplyr::mutate(ID = factor(ID, levels = go_term_sim_order), Count = dplyr::if_else(direction == "up", Count, -Count)) %>%
  dplyr::mutate(celltype = factor(celltype, levels = c("Mono", "CD8T", "CD4T", "NK", "B")))

go_cls_p <- ggplot() +
  geom_segment(aes(x = sqrt(y), y = x, xend = sqrt(yend), yend = xend), data = segment(go_cluster_dd), linewidth = 0.1) +
  geom_rect(aes(ymin = Y.1, ymax = Y.2, xmin = - 0.3 * go_cluster_tree_xmax, xmax = go_cluster_tree_xmax), data = go_cluster_tree_rect, color = "orange", alpha = 0.15, fill = NA) +
  geom_text(aes(x = -0.15 * go_cluster_tree_xmax, y = (Y.1 + Y.2) / 2, label = cluster), data = go_cluster_tree_rect, hjust = 0.5, size = 4, color = rbgc[11]) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_reverse(expand = c(0.02, 0.02)) +
  theme_void()

go_hm_p <- ggplot() +
  geom_tile(aes(x = ID1, y = ID2, fill = similarity), data = go_sim_mat) +
  scale_fill_gradient(low = "#90D5FFEE", high = "darkblue", na.value = "gray99") +
  labs(fill = "Similarity", x = NULL, y = NULL) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  theme_classic() +
  theme(legend.position = "left", axis.text = element_blank(), axis.ticks = element_blank(), axis.line = element_blank())

go_pat_p <- ggplot() +
  geom_tile(aes(x = ID, y = celltype, fill = Count), data = go_pat_tab) +
  facet_grid(direction ~ ., space = "free_y", switch = "y") +
  scale_fill_gradient2(low = rbgc[11], mid = rbgc[6], high = rbgc[1]) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  labs(fill = "Nr. of genes", y = NULL, x = NULL) +
  theme(legend.position = "right", axis.text.x = element_blank(), axis.ticks.x = element_blank(), strip.background = element_blank(), strip.placement = "outside")

go_sim_p <- go_pat_p + go_hm_p + go_cls_p + plot_layout(guides = "collect", design = c(area(1, 1, 1, 6), area(2, 1, 4, 6), area(2, 7, 4, 7)))
go_sim_p_save_to <- file.path(plot_dir, paste0("de_gene.", per_omics, ".all_celltype.BP.go_sim_tree.v2.pdf"))
ggsave(go_sim_p_save_to, go_sim_p, width = 7, height = 6.15)


# GO graph
go_graphml_saveto <- "/home/zzhang/Documents/projects/resources/GeneOntology/go.graphml"
go_graph <- read.graph(go_graphml_saveto, format = "graphml")

selected_go_terms <- go_example_hc %>% dplyr::filter(cluster %in% c("1")) %>% dplyr::pull(label)
vs_idx <- which(as_ids(V(go_graph)) %in% selected_go_terms)
sub_go_graph <- subgraph(go_graph, vs_idx)
shared_selected_node <- c("viral process", "regulation of cell-cell adhesion", "response to interleukin−1", "regulation of cell killing", "cell killing", "integrin−mediated signaling pathway", "cellular defense response") #, "leukocyte proliferation", "regulation of cell-cell adhesion")
monocyte_selected_node <- c("response to tumor necrosis factor", "cytokine−mediated signaling pathway")
sub_go_graph_tbl <- as_tbl_graph(sub_go_graph) %>%
  tidygraph::activate(edges) %>%
  tidygraph::mutate(relationship = ifelse(relationship == "NA" & is_a, "is_a", relationship)) %>%
  tidygraph::activate(nodes) %>%
  # tidygraph::filter(!node_is_isolated()) %>%
  tidygraph::mutate(centrality = centrality_pagerank()) %>%
  # tidygraph::mutate(node_label = ifelse(go_term %in% c(shared_selected_node, monocyte_selected_node), go_term, NA))
  tidygraph::mutate(node_label = ifelse(centrality > 0.025, go_term, ""))
  # %>% tidygraph::mutate(node_label = stringr::str_remove(node_label, "(negative |positive |)regulation of "))

set.seed(31415)
# go_graph_plot <- ggraph(sub_go_graph_tbl, layout = "igraph", algorithm = "stress") +
go_graph_plot <- ggraph(sub_go_graph_tbl, "fr") +
  geom_node_point(aes(size = centrality), alpha = 0.75, color = "black") +
  geom_edge_link(aes(edge_colour = factor(relationship)), start_cap = circle(0.005, 'npc'), end_cap = circle(0.01, 'npc'), arrow = grid::arrow(length = unit(0.01, 'npc'), type="closed")) +
  geom_node_label(aes(label = node_label), segment.linetype = "dotted", fill = "#EEEEEE50", color = "black", size = 5, repel = TRUE, segment.curvature = -0.2, segment.ncp = 2, segment.angle = 30, min.segment.length = 0, max.overlaps = Inf) +
  scale_edge_color_manual(values = c(is_a = "gray90", has_a = "gray90", part_of = "purple", regulates = "orange", negatively_regulates = "blue", positively_regulates = "red")) +
  scale_size(range = c(0.5, 6)) +
  labs(edge_colour = "Relationship", size = "PageRank") +
  theme(legend.position = "right", panel.background = element_blank(), legend.box.background = element_blank())
go_graph_plot_save_to <- file.path(proj_dir, "outputs/analysis/deg/plots", paste0("de_gene.", per_omics, ".all_celltype.BP.go_graph.v2.pdf"))
ggsave(go_graph_plot_save_to, go_graph_plot, width = 8.75, height = 7)


# Number of highly variable genes
citeseq_obj <- readRDS(file.path(proj_dir, "outputs/analysis/integrated/objects/pml_citeseq.rds"))
tar_features <- rownames(citeseq_obj) %>% purrr::discard(~.x %in% c("HBB", "HBA1", "HBA2"))

citeseq_obj <- citeseq_obj[tar_features, ]
DefaultAssay(citeseq_obj) <- "RNA"
citeseq_obj <- FindVariableFeatures(citeseq_obj, selection.method = "vst", nfeatures = 5000)

top_genes <- head(VariableFeatures(citeseq_obj), 40)
p <- VariableFeaturePlot(citeseq_obj, selection.method = "vst")
p <- LabelPoints(plot = p, points = top_genes, repel = TRUE, log = FALSE, max.overlaps = Inf)
p_save_to <- file.path(plot_dir, "highly_variable_genes.citeseq.pdf")
ggsave(p_save_to, plot = p, width=9)


# Functional enrichment of HVGs.
selected_features <- VariableFeatures(citeseq_obj) %>% head(2000)
ego_bp <- enrichGO(selected_features, org.Hs.eg.db, "SYMBOL", "BP")
ego_cc <- enrichGO(selected_features, org.Hs.eg.db, "SYMBOL", "CC")
ego_mf <- enrichGO(selected_features, org.Hs.eg.db, "SYMBOL", "MF")

p_bp <- barplot(ego_bp, showCategory = 12, title="Biological processes", label_format = 40)
p_cc <- barplot(ego_cc, showCategory = 12, title="Cellular components", label_format = 40)
p_mf <- barplot(ego_mf, showCategory = 12, title="Molecular functions", label_format = 40)

p <- p_bp / p_cc / p_mf + plot_layout(guides="auto")
plot_save_to <- file.path(plot_dir, "highly_variable_genes.gene_ontology.citeseq.pdf")
ggsave(plot_save_to, plot = p, height = 15)



# ----------------------------------------------------------------------------------------------------------------------
# similarity
go_sim_p <- ggplot(go_sim_plot_mat) +
  geom_tile(aes(x = ID1, y = ID2, fill = similarity)) +
  scale_fill_gradient(low = rbgc[6], high = rbgc[11]) +
  labs(fill = "GO similarity", x = NULL, y = NULL) +
  theme_classic() +
  theme(legend.position = "left", axis.text = element_blank(), axis.ticks = element_blank(), axis.line = element_blank())

go_p <- go_enr_p + go_sim_p + plot_layout(widths = c(1, 2), guides = "collect") & theme(plot.margin = unit(c(0.01, 0, 0.01, 0), "npc"))
go_plot_save_to <- file.path(plot_dir, paste0("de_gene.", per_omics, ".all_celltype.BP.go_enrichment.pdf"))
ggsave(go_plot_save_to, go_p, width = 8, height = 5.5)


selected_go_tab <- go_tab %>%
  dplyr::filter(ID %in% selected_go_terms, Count >= 3, gene_set %in% c("BL (dw)", "BL (up)")) %>%
  dplyr::group_by(celltype) %>%
  dplyr::slice_min(order_by = p.adjust, n = 1000) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(Description = stringr::str_trunc(Description, 64)) %>%
  dplyr::mutate(Description = forcats::fct_reorder2(Description, log10(p.adjust), celltype))

p <- ggplot(selected_go_tab) +
  geom_point(aes(x = celltype, y = Description, size = Count, color = p.adjust)) +
  scale_color_gradient(high = rbgc[5], low = rbgc[1]) +
  facet_wrap(~gene_set) +
  theme_bw()
selected_go_enr_plot_save_to <- file.path(proj_dir, "outputs/analysis/deg/plots", paste0("de_gene.", per_omics, ".all_celltype.BP.go_enrichment.selected.pdf"))
ggsave(selected_go_enr_plot_save_to, p, width = 7.5, height = 20)



# Dynamics along time
# Dynamic gene list
dym_deg_tab <- de_tab_all %>%
  dplyr::filter(donors == "Responder", p_val_adj < 0.05, abs(avg_log2FC) > 0.25, comparison != "3M.vs.BL", pct.1 > 0.05, pct.2 > 0.05, data_source == "pml_citeseq") %>%
  dplyr::group_by(comparison, celltype)

dym_deg_list <- dym_deg_tab %>%
  dplyr::summarise(gene_set = list(gene_symbol)) %>%
  dplyr::mutate(gene_set_name = paste0(comparison, "(", stringr::str_remove_all(celltype, " "), ")")) %>%
  dplyr::pull(gene_set, gene_set_name)

dym_deg_list %>% lapply(length) %>% unlist() %>% sort

# Interaction of dynamics DE genes
comb_mat <- dym_deg_list %>% list_to_matrix() %>% make_comb_mat() %>% `[`(comb_size(.) >= 5)
gene_upset_save_to <- file.path(proj_dir, "outputs/analysis/deg/plots", paste0("de_gene.", per_omics, ".dynamics.upset_plot.pdf")) %>% stringr::str_remove_all(" ")
pdf(gene_upset_save_to, width = 8, height = 3.75)
UpSet(comb_mat, top_annotation = upset_top_annotation(comb_mat, add_numbers = TRUE), left_annotation = upset_left_annotation(comb_mat, add_numbers = TRUE)) %>% draw()
dev.off()


# Expression pattern of CD8+ T cells
pbmc_int <- readRDS(file.path(object_dir, "pml_citeseq.rds"))
tar_cells <- pbmc_int@meta.data %>%
  dplyr::filter(Response == "Responder", predicted.celltype.l1 == "CD8 T") %>%
  rownames()

expr_pattern_tab <- de_tab_all %>%
  dplyr::filter(donors == "Responder", p_val_adj < 0.05, celltype == "CD8 T", data_source == "pml_citeseq", abs(avg_log2FC) >= 0.25, pct.1 > 0.1, pct.2 > 0.1) %>%
  dplyr::mutate(direction = ifelse(avg_log2FC > 0, "+", "-")) %>%
  dplyr::mutate(comparison = factor(comparison, levels = c("6W.vs.BL", "3M.vs.6W", "3M.vs.BL"))) %>%
  dplyr::group_by(gene_symbol) %>%
  dplyr::arrange(comparison) %>%
  dplyr::summarise(
    comparison_order = paste(comparison, collapse = "|"),
    expression_pattern = dplyr::cur_data() %>% dplyr::mutate(direction = dplyr::if_else(avg_log2FC > 0, "+", "-")) %>% dplyr::pull(direction) %>% paste(collapse = ""),
    expression_pattern = dplyr::case_when(expression_pattern == "---" ~ "Dw", expression_pattern == "+++" ~ "Up", T ~ "Mixed"),
    expression_pattern = factor(expression_pattern, levels = c("Up", "Dw", "Mixed"))
  )

tar_features <- expr_pattern_tab$gene_symbol
norm_avg_expr_tab <- pbmc_int[, tar_cells] %>%
  AverageExpression(features = tar_features, group.by = "Timepoint", assay = "SCT") %>%
  `$`(SCT) %>%
  as.data.frame() %>%
  dplyr::mutate(gene_symbol = rownames(.)) %>%
  dplyr::relocate(gene_symbol)

exp_clusters <- hclust(dist(as.matrix(norm_avg_expr_tab[, -1])), method = "ward.D2")
gene_order <- exp_clusters$labels[exp_clusters$order]

expr_pat_plot_tab <- norm_avg_expr_tab %>%
  dplyr::left_join(expr_pattern_tab, by = "gene_symbol") %>%
  dplyr::select(-c(comparison_order)) %>%
  dplyr::group_by(expression_pattern) %>%
  dplyr::slice_head(n = 20) %>%
  tidyr::pivot_longer(-c(gene_symbol, expression_pattern), names_to = "Timepoint", values_to = "AverageExpression") %>%
  dplyr::mutate(gene_symbol = factor(gene_symbol, levels = gene_order)) %>%
  dplyr::arrange(gene_symbol)

expr_pat_plot_save_to <- file.path(proj_dir, "outputs/analysis/deg", paste0("de_gene.", per_omics, ".dynamics.expression_pattern_heatmap.pdf"))
expr_pat_plot <- ggplot(expr_pat_plot_tab) +
  geom_tile(aes(x = gene_symbol, y = Timepoint, fill = AverageExpression)) +
  scale_fill_gradientn(name = "Avg. exp.", colors = rev(rbgc), limits = c(-10, 10)) +
  # scale_x_discrete(breaks = c("CD2", "ITGB1"), labels = c("CD2", "ITGB1")) +
  labs(x = NULL, y = NULL) +
  facet_grid(. ~ expression_pattern, scales = "free_x", space = "free_x") +
  theme_classic() +
  theme(legend.position = "top", axis.line = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"), panel.spacing = unit(.005, "npc"))
ggsave(expr_pat_plot_save_to, plot = expr_pat_plot, width = 6, height = 3.5)


# Functional enrichment of DEGs at baseline
plot_tab <- go_tab %>%
  dplyr::filter(comparison == "Rs.vs.NRs_BL", data_source == "pml_citeseq", p.adjust <= 0.05, Count >= 3) %>%
  dplyr::group_by(celltype, ONTOLOGY) %>%
  dplyr::slice_min(order_by = pvalue, with_ties = FALSE, n = 5)


# ----------------------------------------------------------
for (pct in tar_cell_types) {
  # GO enrichment of BP, biological processes
  go_plot_tab <- go_tab %>%
    dplyr::filter(ONTOLOGY %in% c("BP"), p.adjust < 0.05, Count >= 5, gene_set %in% c("3M (dw)", "3M (up)", "BL (dw)", "BL (up)")) %>%
    dplyr::group_by(gene_set) %>%
    dplyr::slice_min(order_by = p.adjust, n = 10) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(gene_set = factor(gene_set, levels = c("BL (up)", "3M (up)", "BL (dw)", "3M (dw)"))) %>%
    dplyr::mutate(Short_Description = stringr::str_trunc(Description, 64)) %>%
    dplyr::mutate(Short_Description = forcats::fct_reorder2(Short_Description, log10(p.adjust), gene_set))

  p <- ggplot(go_plot_tab) +
    geom_point(aes(x = gene_set, y = Short_Description, size = -log10(p.adjust), color = log2_odds_ratio)) +
    scale_color_gradient(low = "gray90", high = "red") +
    scale_size(range = c(1, 6)) +
    labs(x = NULL, y = NULL, size = "-log10(p.adj)", color = "log2(OR)") +
    theme_classic() +
    theme(axis.text.y = element_text(size = 12), axis.text.x = element_text(angle = 45, hjust = 1, size = 12))

  go_plot_save_to <- file.path(proj_dir, "outputs/analysis/deg", paste0("de_gene.", per_omics, ".", stringr::str_remove_all(pct, " "), ".go_enrichment.pdf"))
  ggsave(go_plot_save_to, p, width = 8, height = 8)

  # Estimate the similarity of enriched GO terms
  go_sel_tab <- dplyr::filter(go_tab, p.adjust < 0.05, Count >= 5, ONTOLOGY == "BP", gene_set %in% c("3M (dw)", "3M (up)", "BL (dw)", "BL (up)"))
  go_terms <- dplyr::pull(go_sel_tab, ID)
  go_split <- dplyr::pull(go_sel_tab, gene_set)
  go_sim_mat <- mgoSim(go_terms, go_terms, semData = hsGO, measure = "Rel", combine = NULL)
  colnames(go_sim_mat) <- NULL
  rownames(go_sim_mat) <- NULL
  go_sim_mat[is.na(go_sim_mat)] <- 0

  col_fun <- colorRamp2(c(0, max(go_sim_mat)), c("gray95", "red"))
  go_sim_plot_save_to <- file.path(plot_dir, paste0("de_gene.", per_omics, ".", stringr::str_remove_all(pct, " "), ".go_sim_mat.pdf"))
  pdf(go_sim_plot_save_to, width = 8.5, height = 8)
  Heatmap(go_sim_mat, column_split = go_split, row_split = go_split, col = col_fun) %>% draw() # , top_annotation = top_ann)
  dev.off()
}


# UpSet plot to show overlap of DE genes between responders and non-responders at each time point per cell type.
go_tab <- NULL
go_tab_save_to <- file.path(table_dir, paste("de_gene", per_omics, pct, "timepoint_comparison.go_enrichment.csv", sep = ".")) %>% stringr::str_remove_all(" ")
if (file.exists(go_tab_save_to)) {
  cat(paste("Reading", go_tab_save_to, "...\n"))
  go_tab <- fread(go_tab_save_to) %>% rbind(go_tab)
} else {
  for (pct in tar_cell_types) {
    # pct <- "CD8 T"
    selected_tab <- de_tab_all %>%
      dplyr::filter(celltype == pct, data_source == per_omics, donors == "Responder", p_val_adj < 0.05, abs(avg_log2FC) > 0.25) %>%
      dplyr::select(gene_symbol, direction, comparison)

    if (nrow(selected_tab) == 0) { next }

    tryCatch({
      gene_set_list <- selected_tab %>%
        (function(tab) {
          list(
            `6W.vs.BL (up)` = tab %>% dplyr::filter(direction == "up", comparison == "6W.vs.BL") %>% dplyr::pull(gene_symbol),
            `3M.vs.6W (up)` = tab %>% dplyr::filter(direction == "up", comparison == "3M.vs.6W") %>% dplyr::pull(gene_symbol),
            `3M.vs.BL (up)` = tab %>% dplyr::filter(direction == "up", comparison == "3M.vs.BL") %>% dplyr::pull(gene_symbol),
            `6W.vs.BL (dw)` = tab %>% dplyr::filter(direction == "dw", comparison == "6W.vs.BL") %>% dplyr::pull(gene_symbol),
            `3M.vs.6W (dw)` = tab %>% dplyr::filter(direction == "dw", comparison == "3M.vs.6W") %>% dplyr::pull(gene_symbol),
            `3M.vs.BL (dw)` = tab %>% dplyr::filter(direction == "dw", comparison == "3M.vs.BL") %>% dplyr::pull(gene_symbol)
        )})

      # GO enrichment
      per_go_tab <- lapply(names(gene_set_list), function(x, .gsl) {
        gene_vec <- .gsl[[x]]
        tmp <- enrichGO(gene_vec, org.Hs.eg.db, "SYMBOL", "ALL")
        tmp@result %>% dplyr::mutate(gene_set = x, label = paste0(x, " (n=", length(gene_vec), ")"))
      }, .gsl = gene_set_list) %>%
        Reduce(rbind, .) %>%
        dplyr::mutate(bg_ratio_dec = lapply(.$BgRatio, function(x) eval(parse(text = x))) %>% unlist) %>%
        dplyr::mutate(gene_ratio_dec = lapply(.$GeneRatio, function(x) eval(parse(text = x))) %>% unlist) %>%
        dplyr::mutate(log2_odds_ratio = log2(gene_ratio_dec / bg_ratio_dec))

      go_tab <- rbind(go_tab, per_go_tab)

      # # Number of identified genes
      # gene_upset_save_to <- file.path(plot_dir, paste("de_gene", per_omics, pct, "timepoint_comparison.upset_plot.pdf", sep = ".")) %>% stringr::str_remove_all(" ")
      # if (file.exists(gene_upset_save_to)) {
      #   cat("File", gene_upset_save_to, "exists. Skipping ...\n")
      # } else {
      #   comb_mat <- gene_set_list %>% list_to_matrix() %>% make_comb_mat()
      #   pdf(gene_upset_save_to, width = 6, height = 3.25)
      #   UpSet(comb_mat) %>% draw()
      #   dev.off()
      # }

      # # GO plot
      # go_plot_save_to <- file.path(plot_dir, paste("de_gene", per_omics, pct, "timepoint_comparison.go_enrichment.pdf", sep = ".")) %>% stringr::str_remove_all(" ")
      # if (file.exists(go_plot_save_to)) {
      #   cat("GO plot exists\n")
      # } else {
      #   go_plot_tab <- per_go_tab %>%
      #     dplyr::filter(ONTOLOGY %in% c("BP"), p.adjust < 0.05, Count >= 5, gene_set %in% names(gene_set_list)) %>%
      #     dplyr::group_by(gene_set) %>%
      #     dplyr::slice_min(order_by = p.adjust, n = 10) %>%
      #     dplyr::ungroup() %>%
      #     dplyr::mutate(gene_set = factor(gene_set, levels = names(gene_set_list))) %>%
      #     dplyr::mutate(Short_Description = stringr::str_trunc(Description, 64)) %>%
      #     dplyr::mutate(Short_Description = forcats::fct_reorder2(Short_Description, log10(p.adjust), gene_set))

      #   p <- ggplot(go_plot_tab) +
      #     geom_point(aes(x = gene_set, y = Short_Description, size = -log10(p.adjust), color = log2_odds_ratio)) +
      #     scale_color_gradient(low = "gray90", high = "red") +
      #     scale_size(range = c(1, 6)) +
      #     labs(x = NULL, y = NULL, size = "-log10(p.adj)", color = "log2(OR)") +
      #     theme_classic() +
      #     theme(axis.text.y = element_text(size = 12), axis.text.x = element_text(angle = 45, hjust = 1, size = 12))
      #   ggsave(go_plot_save_to, p, width = 8, height = 8)
      # }

      # # GO terms similarity 
      # go_sim_plot_save_to <- file.path(plot_dir, paste("de_gene", per_omics, pct, "timepoint_comparison.go_sim_mat.pdf", sep = ".")) %>% stringr::str_remove_all(" ")
      # if (file.exists(go_sim_plot_save_to)) {
      #   cat("GO similarity plot exists\n")
      # } else {
      #   go_sel_tab <- dplyr::filter(per_go_tab, p.adjust < 0.05, Count >= 5, ONTOLOGY == "BP") %>% dplyr::group_by(gene_set) %>% dplyr::slice_min(n = 30, order_by = p.adjust)
      #   go_terms <- dplyr::pull(go_sel_tab, ID)
      #   go_split <- dplyr::pull(go_sel_tab, gene_set)
      #   go_sim_mat <- mgoSim(go_terms, go_terms, semData = hsGO, measure = "Rel", combine = NULL)
      #   colnames(go_sim_mat) <- NULL
      #   rownames(go_sim_mat) <- NULL
      #   go_sim_mat[is.na(go_sim_mat)] <- 0

      #   col_fun <- colorRamp2(c(0, max(go_sim_mat)), c("gray95", "red"))
      #   pdf(go_sim_plot_save_to, width = 8.5, height = 8)
      #   Heatmap(go_sim_mat, column_split = go_split, row_split = go_split, col = col_fun)# , top_annotation = top_ann)
      #   dev.off()
      # }
    }, error = function(e) cat(e$message, paste(pct, "failed\n")))
  }

  fwrite(go_tab, go_tab_save_to)
}
