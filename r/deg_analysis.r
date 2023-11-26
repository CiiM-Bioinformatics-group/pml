#!/usr/bin/env Rscript
# File: deg_analysis.r
# Author: Zhenhua Zhang
# E-mail: zhenhua.zhang217@gmail.com
# Created: Jun 28, 2023
# Updated: Nov 23, 2023

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(Seurat)

  library(ggsci)
  library(ggbreak)
  library(ggrepel)
  library(ggsankey) # or library(ggalluvial)
  library(ggVennDiagram)

  library(gtools)

  library(GOplot)
  library(enrichplot)
  library(ComplexHeatmap)

  library(GOSemSim)

  library(circlize)

  library(clusterProfiler)
  library(org.Hs.eg.db)
})


#
## DEG identification
#
proj_dir <- "~/Documents/projects/wp_pml"
deg_dir <- file.path(proj_dir, "outputs/analysis/deg")
plot_dir <- file.path(proj_dir, "outputs/analysis/integrated/plots")
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

de_tab_save_to <- file.path(proj_dir, "outputs/analysis/deg", paste0("pbmc.integrated.de_gene.csv"))
cp_tab_save_to <- file.path(proj_dir, "outputs/analysis/deg", paste0("pbmc.integrated.cell_proportion.csv"))
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
        dplyr::select(PatientID = Vireo_assignment, Timepoint, Response, predicted.celltype.l1, prop)

      ## 1. Cell proportion overview
      p <- ggplot(cell_proportion_tab) +
        geom_line(aes(x = Timepoint, y = prop, group = PatientID)) +
        geom_point(aes(x = Timepoint, y = prop, color = Response), alpha = 0.75) +
        scale_fill_npg() +
        theme_classic() +
        theme(axis.text.y = element_text(size = 8), axis.title = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(y = "Proportion", x = "Cell type") +
        facet_grid(~predicted.celltype.l1, space = "free_x")
      save_to <- file.path(proj_dir, "outputs/analysis/cell_proportion", paste0("pbmc.", omics_type, ".integrated.cell_proportion.pdf"))
      ggsave(save_to, plot = p, width = 8.5, height = 3)

      ## 2. Cell proportion by patient
      p <- ggplot(cell_proportion_tab) +
        geom_point(aes(x = PatientID, y = prop, color = predicted.celltype.l1)) +
        scale_color_npg() +
        theme_classic() +
        theme(axis.text.y = element_text(size = 8), axis.title = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(y = "Proportion", x = "Cell type") +
        facet_grid(~Timepoint, scales = "free_x", space = "free_x")
      save_to <- file.path(proj_dir, "outputs/analysis/cell_proportion", paste0("pbmc.", omics_type, ".integrated.cell_proportion_by_patient.pdf"))
      ggsave(save_to, plot = p, width = 8.5, height = 3)

      ## 3. Cell proportion by timepoints
      p <- ggplot(cell_proportion_tab) +
        geom_col(aes(x = Timepoint, y = prop, fill = predicted.celltype.l1), width = 0.985) +
        scale_fill_npg() +
        theme_classic() +
        theme(axis.text.y = element_text(size = 8), axis.title = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(y = "Proportion", x = "Cell type") +
        facet_grid(~PatientID, scales = "free_x", space = "free_x")
      save_to <- file.path(proj_dir, "outputs/analysis/cell_proportion", paste0("pbmc.", omics_type, ".integrated.cell_proportion_by_timepoints.pdf"))
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
      save_to <- file.path(proj_dir, "outputs/analysis/cell_proportion", paste0("pbmc.", omics_type, ".integrated.cell_proportion.pdf"))
      ggsave(save_to, plot = p, width = 8.5, height = 3)

      # 2. Cell proportion by patient
      p <- ggplot(cell_proportion_tab) +
        geom_point(aes(x = PatientID, y = prop, color = predicted.celltype.l1)) +
        scale_color_npg() +
        theme_classic() +
        theme(axis.text.y = element_text(size = 8), axis.title = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(y = "Proportion", x = "Cell type") +
        facet_grid(~Timepoint, scales = "free_x", space = "free_x")
      save_to <- file.path(proj_dir, "outputs/analysis/cell_proportion", paste0("pbmc.", omics_type, ".integrated.cell_proportion_by_patient.pdf"))
      ggsave(save_to, plot = p, width = 8.5, height = 3)

      # 3. Cell proportion by timepoints
      p <- ggplot(cell_proportion_tab) +
        geom_col(aes(x = Timepoint, y = prop, fill = predicted.celltype.l1), width = 0.985) +
        scale_fill_npg() +
        theme_classic() +
        theme(axis.text.y = element_text(size = 8), axis.title = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(y = "Proportion", x = "Cell type") +
        facet_grid(~PatientID, scales = "free_x", space = "free_x")
      save_to <- file.path(proj_dir, "outputs/analysis/cell_proportion", paste0("pbmc.", omics_type, ".integrated.cell_proportion_by_timepoints.pdf"))
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
    save_to <- file.path(proj_dir, "outputs/analysis/cell_proportion", paste0("pbmc.", omics_type, "integrated.cell_propotion.csv"))
    fwrite(cell_proportion_tab, save_to)

    de_tab <- dplyr::mutate(de_tab, direction = dplyr::if_else(avg_log2FC > 0, "up", "dw"), data_source = omics_type)
    save_to <- file.path(proj_dir, "outputs/analysis/deg", paste0("pbmc.", omics_type, ".integrated.de_gene.csv"))
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
per_omics <- "pml_citeseq"
de_fmt_tab <- as.data.table(de_tab_all) %>%
  dplyr::filter(p_val_adj <= 0.05, donors == "All", data_source == per_omics, abs(avg_log2FC) >= 0.25) %>% 
  tidyr::separate(comparison, into = c("comparison", "timepoint"), sep = "_") %>%
  dplyr::mutate(timepoint = factor(timepoint, levels = c("BL", "3M", "PT")))

# Bar plot show number of DE genes
de_count_tab <- de_fmt_tab %>%
  dplyr::filter(data_source == per_omics) %>%
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
  facet_wrap(~timepoint, ncol = 1, strip.position = "right")

save_to <- file.path(proj_dir, "outputs/analysis/deg", paste0("de_gene.", per_omics, ".Rs_vs_NR.per_timepoint.pdf"))
ggsave(save_to, plot = p, width = 4, height = 12)

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
  save_to <- file.path(proj_dir, "outputs/analysis/deg", paste("de_gene", per_omics, pct, "Rs_vs_NR.per_timepoint.venn_diagram.v2.pdf", sep = "."))
  p <- ggVennDiagram(venn_plot_list) + scale_fill_gradient(low = "white", high = "darkblue") + theme(legend.position = "none")
  ggsave(save_to, plot = p, width = 6, height = 6)
}


# Up regulated genes in monocytes and CD8+ T, two time points
pct <- "CD4 T"
pct <- "CD8 T"
pct <- "Mono"

for (pct in tar_cell_types) {
  save_to <- file.path(proj_dir, "outputs/analysis/deg", paste0("de_gene.", per_omics, ".", stringr::str_remove_all(pct, " "), ".Rs_vs_NR.go_enrichment.csv"))
  if (!file.exists(save_to)) {
    comb_mat <- de_fmt_tab %>%
      dplyr::filter(celltype == pct, data_source == per_omics, p_val_adj < 0.05, abs(avg_log2FC) > 0.25) %>%
      dplyr::select(gene_symbol, direction, timepoint) %>%
      tidyr::pivot_wider(values_from = direction, names_from = timepoint) %>%
      tidyr::pivot_longer(c(BL, `3M`)) %>%
      dplyr::filter(!is.na(value)) %>%
      (function(tab) {
        list(
          `BL (up)` = tab %>% dplyr::filter(value == "up", name == "BL") %>% dplyr::pull(gene_symbol),
          `3M (up)` = tab %>% dplyr::filter(value == "up", name == "3M") %>% dplyr::pull(gene_symbol),
          `BL (dw)` = tab %>% dplyr::filter(value == "dw", name == "BL") %>% dplyr::pull(gene_symbol),
          `3M (dw)` = tab %>% dplyr::filter(value == "dw", name == "3M") %>% dplyr::pull(gene_symbol)
        )}) %>%
      list_to_matrix() %>%
      make_comb_mat()

    gene_set_list <- list(
      `BL (up)` = lapply(c("1100", "1001", "1000"), extract_comb, m = comb_mat) %>% unlist(),
      `3M (up)` = lapply(c("1100", "0110", "0100"), extract_comb, m = comb_mat) %>% unlist(),
      `BL (dw)` = lapply(c("0110", "0011", "0010"), extract_comb, m = comb_mat) %>% unlist(),
      `3M (dw)` = lapply(c("1001", "0111", "0001"), extract_comb, m = comb_mat) %>% unlist(),
      `Shared` = lapply(c("1100"), extract_comb, m = comb_mat) %>% unlist(),
      `3M (only)` = lapply(c("0100", "0110"), extract_comb, m = comb_mat) %>% unlist(),
      `BL (only)` = lapply(c("1000", "1001"), extract_comb, m = comb_mat) %>% unlist()
    )

    go_tab <- lapply(names(gene_set_list), function(x, .gsl) {
      gene_vec <- .gsl[[x]]
      tmp <- enrichGO(
        gene = gene_vec, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "ALL", pAdjustMethod = "BH", pvalueCutoff = 0.05,
        qvalueCutoff = 0.05, readable = TRUE
      )
      tmp@result %>% dplyr::mutate(gene_set = x, label = paste0(x, " (n=", length(gene_vec), ")"))
    }, .gsl = gene_set_list) %>%
      Reduce(rbind, .) %>%
      dplyr::mutate(bg_ratio_dec = lapply(.$BgRatio, function(x) eval(parse(text = x))) %>% unlist) %>%
      dplyr::mutate(gene_ratio_dec = lapply(.$GeneRatio, function(x) eval(parse(text = x))) %>% unlist) %>%
      dplyr::mutate(log2_odds_ratio = log2(gene_ratio_dec / bg_ratio_dec))

    fwrite(go_tab, save_to)
  } else {
    go_tab <- fread(save_to)
  }

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
  go_sim_plot_save_to <- file.path(proj_dir, "outputs/analysis/deg", paste0("de_gene.", per_omics, ".", stringr::str_remove_all(pct, " "), ".go_sim_mat.pdf"))
  pdf(go_sim_plot_save_to, width = 8.5, height = 8)
  Heatmap(go_sim_mat, column_split = go_split, row_split = go_split, col = col_fun)# , top_annotation = top_ann)
  dev.off()

  # # GO similarity score
  # go_simscore_mat <- go_sel_tab$gene_set %>% unique %>% permutations(4, 2, v = ., repeats.allowed = TRUE) %>% as.data.frame() %>% apply(1, function(x, tab) {
  #   id_set_1 <- tab %>% dplyr::filter(gene_set == x[1]) %>% dplyr::pull(ID)
  #   id_set_2 <- tab %>% dplyr::filter(gene_set == x[2]) %>% dplyr::pull(ID)
  #   go_sim_score <- mgoSim(id_set_1, id_set_2, semData = hsGO, measure = "Jiang", combine = "BMA")
  #   data.frame(set_1 = x[1], set_2 = x[2], GO_SimScore = go_sim_score)
  # }, tab = go_sel_tab) %>%
  #   Reduce(rbind, .) %>%
  #   tidyr::pivot_wider(names_from = set_2, values_from = GO_SimScore) %>%
  #   tibble::column_to_rownames("set_1") %>%
  #   as.matrix()
  # save_to <- file.path(proj_dir, "outputs/analysis/deg", paste0("de_gene.", per_omics, ".", stringr::str_remove_all(pct, " "), ".go_sim_score.pdf"))
  # pdf(save_to, width = 4.5, height = 4)
  # Heatmap(go_simscore_mat, col = col_fun)# , top_annotation = top_ann)
  # dev.off()
}


# UpSet plot to show overlap of DE genes between responders and non-responders at each time point per cell type.
for (pct in tar_cell_types) {
  # pct <- "CD4 T"
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

    # Number of identified genes
    comb_mat <- gene_set_list %>% list_to_matrix() %>% make_comb_mat()
    gene_upset_save_to <- file.path(proj_dir, "outputs/analysis/deg", paste("de_gene", per_omics, pct, "timepoint_comparison.upset_plot.pdf", sep = ".")) %>% stringr::str_remove_all(" ")
    pdf(gene_upset_save_to, width = 6, height = 3.25)
    UpSet(comb_mat) %>% draw()
    dev.off()

    # GO enrichment
    go_tab_save_to <- file.path(proj_dir, "outputs/analysis/deg", paste("de_gene", per_omics, pct, "timepoint_comparison.go_enrichment.csv", sep = ".")) %>% stringr::str_remove_all(" ")
    if (file.exists(go_tab_save_to)) {
      cat(paste("Reading", go_tab_save_to, "...\n"))
      go_tab <- fread(go_tab_save_to)
    } else {
      cat(paste("Calculating", go_tab_save_to, "...\n"))
      go_tab <- lapply(names(gene_set_list), function(x, .gsl) {
        gene_vec <- .gsl[[x]]
        tmp <- enrichGO(
          gene = gene_vec, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "ALL", pAdjustMethod = "BH", pvalueCutoff = 0.05,
          qvalueCutoff = 0.05, readable = TRUE
        )
        tmp@result %>% dplyr::mutate(gene_set = x, label = paste0(x, " (n=", length(gene_vec), ")"))
      }, .gsl = gene_set_list) %>%
        Reduce(rbind, .) %>%
        dplyr::mutate(bg_ratio_dec = lapply(.$BgRatio, function(x) eval(parse(text = x))) %>% unlist) %>%
        dplyr::mutate(gene_ratio_dec = lapply(.$GeneRatio, function(x) eval(parse(text = x))) %>% unlist) %>%
        dplyr::mutate(log2_odds_ratio = log2(gene_ratio_dec / bg_ratio_dec))
      fwrite(go_tab, go_tab_save_to)
    }

    # GO plot
    go_plot_save_to <- file.path(proj_dir, "outputs/analysis/deg", paste("de_gene", per_omics, pct, "timepoint_comparison.go_enrichment.pdf", sep = ".")) %>% stringr::str_remove_all(" ")
    if (file.exists(go_plot_save_to)) {
      cat("GO plot exists\n")
    } else {
      go_plot_tab <- go_tab %>%
        dplyr::filter(ONTOLOGY %in% c("BP"), p.adjust < 0.05, Count >= 5, gene_set %in% names(gene_set_list)) %>%
        dplyr::group_by(gene_set) %>%
        dplyr::slice_min(order_by = p.adjust, n = 10) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(gene_set = factor(gene_set, levels = names(gene_set_list))) %>%
        dplyr::mutate(Short_Description = stringr::str_trunc(Description, 64)) %>%
        dplyr::mutate(Short_Description = forcats::fct_reorder2(Short_Description, log10(p.adjust), gene_set))
      p <- ggplot(go_plot_tab) +
        geom_point(aes(x = gene_set, y = Short_Description, size = -log10(p.adjust), color = log2_odds_ratio)) +
        scale_color_gradient(low = "gray90", high = "red") +
        scale_size(range = c(1, 6)) +
        labs(x = NULL, y = NULL, size = "-log10(p.adj)", color = "log2(OR)") +
        theme_classic() +
        theme(axis.text.y = element_text(size = 12), axis.text.x = element_text(angle = 45, hjust = 1, size = 12))
      ggsave(go_plot_save_to, p, width = 8, height = 8)
    }

    # GO terms similarity 
    go_sim_plot_save_to <- file.path(proj_dir, "outputs/analysis/deg", paste("de_gene", per_omics, pct, "timepoint_comparison.go_sim_mat.pdf", sep = ".")) %>% stringr::str_remove_all(" ")
    if (file.exists(go_sim_plot_save_to)) {
      cat("GO similarity plot exists\n")
    } else {
      go_sel_tab <- dplyr::filter(go_tab, p.adjust < 0.05, Count >= 5, ONTOLOGY == "BP") %>% dplyr::group_by(gene_set) %>% dplyr::slice_min(n = 30, order_by = p.adjust)
      go_terms <- dplyr::pull(go_sel_tab, ID)
      go_split <- dplyr::pull(go_sel_tab, gene_set)
      go_sim_mat <- mgoSim(go_terms, go_terms, semData = hsGO, measure = "Rel", combine = NULL)
      colnames(go_sim_mat) <- NULL
      rownames(go_sim_mat) <- NULL
      go_sim_mat[is.na(go_sim_mat)] <- 0

      col_fun <- colorRamp2(c(0, max(go_sim_mat)), c("gray95", "red"))
      pdf(go_sim_plot_save_to, width = 8.5, height = 8)
      Heatmap(go_sim_mat, column_split = go_split, row_split = go_split, col = col_fun)# , top_annotation = top_ann)
      dev.off()
    }

  }, error = function(e) cat(e$message, paste(pct, "failed\n")) )
}


# DE genes replicated by RNA-seq
dis_deg_tab <- de_tab_all %>%
  dplyr::filter(donors == "All", comparison == "Rs.vs.NRs_BL", data_source == "pml_citeseq", p_val_adj < 0.05, abs(avg_log2FC) > 0.25) %>%
  dplyr::select(celltype, gene_symbol, direction, p_val, avg_log2FC)

val_deg_tab <- de_tab_all %>%
  dplyr::filter(donors == "All", comparison == "Rs.vs.NRs_BL", data_source == "pml_rnaseq", p_val < 0.05) %>%
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
save_to <- file.path(proj_dir, "outputs/analysis/deg/de_gene.replicated_at_BL.pdf")
ggsave(save_to, plot = p, width = 4, height = 4.5)






# sessionInfo()
# R version 4.2.0 (2022-04-22)
# Platform: x86_64-pc-linux-gnu (64-bit)                          
# Running under: AlmaLinux 8.5 (Arctic Sphynx)                     
#                                                                   
# Matrix products: default                                                
# BLAS:   /vol/projects/zzhang/tools/R/4.2.0/lib64/R/lib/libRblas.so
# LAPACK: /vol/projects/zzhang/tools/R/4.2.0/lib64/R/lib/libRlapack.so    
#                                                                     
# locale:                                                         
# [1] C                                                           
#                                                                    
# attached base packages:                                          
# [1] stats4    grid      stats     graphics  grDevices utils     datasets
# [8] methods   base                                             
#                                                                 
# other attached packages:                                           
#  [1] org.Hs.eg.db_3.15.0   AnnotationDbi_1.58.0  IRanges_2.30.0   
#  [4] S4Vectors_0.34.0      Biobase_2.56.0        BiocGenerics_0.42.0
#  [7] clusterProfiler_4.4.4 circlize_0.4.15       GOSemSim_2.22.0
# [10] ComplexHeatmap_2.12.1 enrichplot_1.16.1     GOplot_1.0.2    
# [13] RColorBrewer_1.1-3    gridExtra_2.3         ggdendro_0.1.23
# [16] gtools_3.9.4          ggVennDiagram_1.2.2   ggsankey_0.0.99999
# [19] ggrepel_0.9.1         ggbreak_0.1.1         ggsci_2.9       
# [22] sp_1.5-0              SeuratObject_4.1.0    Seurat_4.1.1        
# [25] data.table_1.14.2     forcats_0.5.1         stringr_1.4.1
# [28] dplyr_1.0.9           purrr_0.3.4           readr_2.1.2    
# [31] tidyr_1.2.0           tibble_3.1.8          ggplot2_3.4.0    
# [34] tidyverse_1.3.2                                            
#                                                                
# loaded via a namespace (and not attached):                     
#   [1] utf8_1.2.2             reticulate_1.25        tidyselect_1.2.0
#   [4] RSQLite_2.3.1          htmlwidgets_1.5.4      BiocParallel_1.30.3
#   [7] Rtsne_0.16             scatterpie_0.1.7       munsell_0.5.0
#  [10] codetools_0.2-18       ica_1.0-3              future_1.26.1        
#  [13] miniUI_0.1.1.1         withr_2.5.0            spatstat.random_3.0-1
#  [16] colorspace_2.0-3       progressr_0.10.1       ROCR_1.0-11      
#  [19] tensor_1.5             DOSE_3.22.0            listenv_0.8.0   
#  [22] GenomeInfoDbData_1.2.8 polyclip_1.10-0        farver_2.1.1   
#  [25] bit64_4.0.5            downloader_0.4         treeio_1.20.1 
#  [28] parallelly_1.32.0      vctrs_0.5.1            generics_0.1.3
#  [31] doParallel_1.0.17      R6_2.5.1               GenomeInfoDb_1.35.15
#  [34] clue_0.3-61            graphlayouts_0.8.0     RVenn_1.1.0 
#  [37] fgsea_1.22.0           bitops_1.0-7           spatstat.utils_3.0-1
#  [40] cachem_1.0.6           gridGraphics_0.5-1     assertthat_0.2.1
#  [43] promises_1.2.0.1       scales_1.2.1           ggraph_2.0.5
#  [46] googlesheets4_1.0.0    rgeos_0.5-9            gtable_0.3.0 
#  [49] globals_0.15.1         goftest_1.2-3          tidygraph_1.2.1
#  [52] rlang_1.0.6            GlobalOptions_0.1.2    splines_4.2.0
#  [55] lazyeval_0.2.2         gargle_1.2.0           spatstat.geom_3.0-3
#  [58] broom_1.0.0            reshape2_1.4.4         abind_1.4-5
#  [61] modelr_0.1.8           backports_1.4.1        httpuv_1.6.5
#  [64] qvalue_2.15.0          tools_4.2.0            ggplotify_0.1.0
#  [67] ellipsis_0.3.2         spatstat.core_2.4-4    ggridges_0.5.3
#  [70] Rcpp_1.0.10            plyr_1.8.7             zlibbioc_1.42.0
#  [73] RCurl_1.98-1.10        rpart_4.1.16           deldir_1.0-6
#  [76] GetoptLong_1.0.5       viridis_0.6.2          pbapply_1.5-0
#  [79] cowplot_1.1.1          zoo_1.8-10             haven_2.5.0
#  [82] cluster_2.1.3          fs_1.5.2               magrittr_2.0.3
#  [85] scattermore_0.8        DO.db_2.9              lmtest_0.9-40
#  [88] reprex_2.0.1           RANN_2.6.1             googledrive_2.0.0
#  [91] fitdistrplus_1.1-8     matrixStats_0.62.0     hms_1.1.1
#  [94] patchwork_1.1.1        mime_0.12              xtable_1.8-4
#  [97] readxl_1.4.0           shape_1.4.6            compiler_4.2.0
# [100] shadowtext_0.1.2       KernSmooth_2.23-20     crayon_1.5.1
# [103] htmltools_0.5.3        ggfun_0.0.6            mgcv_1.8-40
# [106] later_1.3.0            tzdb_0.3.0             aplot_0.1.6
# [109] lubridate_1.8.0        DBI_1.1.3              tweenr_1.0.2
# [112] dbplyr_2.2.1           MASS_7.3-56            Matrix_1.4-1
# [115] cli_3.4.1              parallel_4.2.0         igraph_1.3.4
# [118] pkgconfig_2.0.3        plotly_4.10.0          spatstat.sparse_3.0-0
# [121] foreach_1.5.2          xml2_1.3.3             ggtree_3.4.1
# [124] XVector_0.36.0         rvest_1.0.2            yulab.utils_0.0.5
# [127] digest_0.6.29          sctransform_0.3.3      RcppAnnoy_0.0.19
# [130] spatstat.data_3.0-0    Biostrings_2.64.0      fastmatch_1.1-3
# [133] cellranger_1.1.0       leiden_0.4.2           tidytree_0.3.9
# [136] uwot_0.1.14            shiny_1.7.2            rjson_0.2.21
# [139] lifecycle_1.0.3        nlme_3.1-157           jsonlite_1.8.0
# [142] viridisLite_0.4.0      fansi_1.0.3            pillar_1.8.1
# [145] lattice_0.20-45        KEGGREST_1.36.3        fastmap_1.1.0
# [148] httr_1.4.3             survival_3.3-1         GO.db_3.15.0
# [151] glue_1.6.2             iterators_1.0.14       png_0.1-7
# [154] bit_4.0.4              ggforce_0.3.3          stringi_1.7.8
# [157] blob_1.2.3             memoise_2.0.1          ape_5.6-2
# [160] irlba_2.3.5            future.apply_1.9.0
