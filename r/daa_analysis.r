#!/usr/bin/env Rscript
# File: daa_analysis.r
# Author: Zhenhua Zhang
# E-mail: zhenhua.zhang217@gmail.com
# Created: Nov 02, 2023
# Updated: Nov 26, 2023

suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(patchwork)
  library(dichromat)
  library(ggsci)
  library(ggsignif)
  library(ggrepel)

  library(org.Hs.eg.db)
  library(enrichplot)
  library(RColorBrewer)
  library(clusterProfiler)
  library(igraph)

  # GSEA
  library(fgsea)
  library(msigdbr)

  library(miloR)
  library(SingleCellExperiment)
  library(Seurat)
  library(scran)
  library(scater)
})


rm(list = ls()); gc()

rbgc <- brewer.pal(n = 11, name = "RdBu")
npg3 <- pal_npg("nrc")(3)
lacet3 <- pal_lancet("lanonc")(3)

RANDOM_SEED <- 31415926
overwrite <- TRUE
proj_dir <- "~/Documents/projects/wp_pml"
plot_dir <- file.path(proj_dir, "outputs/analysis/differential_abundance/plots")
object_dir <- file.path(proj_dir, "outputs/analysis/differential_abundance/objects")

my_colors <- colorRampPalette(c("red"))(10)
chrom_xy_genes <- file.path(proj_dir, "inputs/references/sextual_chromosome_genes.txt") %>% fread(header = F) %>% pull(V2)

important_markers <- c(
  c("PDCD1", "CTLA4", "LAG3", "HAVCR2", "CD244", "CD160", "TIGIT"), # T cell exaustion markers, highly expressed inhibitory receptors. HAVCR2 aka TIM-3, CD244 aka 2B4/SLAMF4
  c("TNF", "IL2"), # Lower expression upon T cell exaustion
  c("KLRG1", "IL7R", "CX3CR1", "TBX21", "BACH2"), # High cytoxcity, high proliferation CD8+ Tem, PMID: 29625895
  c("CCR7", "SELL", "CD27", "CD28"), # Maker of T-cell activation in response to antigion stimulation, PMID: 17371966
  # c("TCF7", "LEF1"), # Constant supervision to mature CD8+ T cells
  c("LAMP1"), # Markers for degranulation in activated CD8+ T cells.

  # Cytotoxic CD8+ T lymphocytes mediate the killing of target cells via two major pathways:
  c("PRF1", "GZMA", "GZMB", "GZMH", "GNLY"), # (1) Perforin-granzyme-mediated activation of programmed cell death or apoptosis
  c("FASLG") # (2) fas–fas ligand–mediated induction of apoptosis.
)


# per_omics <- "pml_rnaseq"; pca_source <- "PCA"; umap_source <- "UMAP"
# per_omics <- "pml_rnaseq"; pca_source <- "REF.SPCA"; umap_source <- "REF.UMAP"
# per_omics <- "pml_citeseq"; pca_source <- "REF.SPCA"; umap_source <- "REF.UMAP"
per_omics <- "pml_citeseq"; pca_source <- "INTPCA"; umap_source <- "WNN.UMAP"
pca_source_dict <- c("PCA" = "internal_pca", "INTPCA" = "internal_pca", "REF.SPCA" = "external_pca")

pbmc_int <- file.path(proj_dir, "outputs/analysis/integrated/objects", paste0(per_omics, ".rds")) %>% readRDS()
pbmc_int_sce <- as.SingleCellExperiment(pbmc_int, assay = "SCT")
pbmc_int_sce$Response <- ifelse(pbmc_int_sce$Response == "Non-responder", "Nonresponder", "Responder") %>% factor(levels = c("Nonresponder", "Responder"))
pbmc_int_sce$Sample_by_Timepoint <- paste0(pbmc_int_sce$Vireo_assignment, "_", pbmc_int_sce$Timepoint)
pbmc_int_sce$Treatment <- ifelse(pbmc_int_sce$Timepoint == "BL", "BL", "PT") %>% factor(levels = c("BL", "PT"))

# Check the cell type clusters using UMAPs 
if (FALSE) {
  save_to <- file.path(plot_dir, paste("daa", per_omics, "ref_umap.scran.pdf", sep = "."))
  p <- plotReducedDim(pbmc_int_sce, "REF.UMAP", colour_by = "Response")
  ggsave(save_to, plot = p, width = 5, height = 5)

  save_to <- file.path(plot_dir, paste("daa", per_omics, "umap.scran.pdf", sep = "."))
  p <- plotReducedDim(pbmc_int_sce, umap_source, colour_by = "Response")
  ggsave(save_to, plot = p, width = 5, height = 5)
}


#
## Differential Abundance Analysis
#
alpha <- 0.1 # Significance level used for SpatialFDR
reductdim_d <- 30 # Number of demensions (e.g., PC) used in the neighborhood
nhood_cell_prop <- 0.5 # Proportion of cells used to build the neighborhood
nhood_k <- 30 # Number of neighbors for a cell.
# for (per_comp in c("Responder_by_time")) {
for (per_comp in c("Responder_by_time", "Rs_vs_NRs")) {
  if (per_comp == "Rs_vs_NRs") {
    tar_cells <- pbmc_int_sce %>% colData() %>% as.data.frame() %>% dplyr::filter(Timepoint == "BL") %>% rownames()
    umap_color_by <- "Response"
  } else if (per_comp == "Responder_by_time") {
    tar_cells <- pbmc_int_sce %>% colData() %>% as.data.frame() %>% dplyr::filter(Response == "Responder", predicted.celltype.l1 %in% c("CD8 T", "Mono")) %>% rownames()
    umap_color_by <- "Timepoint"
  }

  save_token <- paste("daa", per_omics, per_comp, pca_source_dict[pca_source], "v2", sep = ".") # Token to distinguish different results
  obj_save_to <- file.path(object_dir, paste0(save_token, ".da_test_object.rds"))
  da_result_save_to <- file.path(object_dir, paste0(save_token, ".da_test_results.csv"))
  da_deg_result_save_to <- file.path(object_dir, paste0(save_token, ".da_deg_test_results.csv"))
  if (all(file.exists(c(obj_save_to, da_result_save_to, da_deg_result_save_to)))) {
    cat("[I]: Loading saved objects ...\n")
    pml_milo <- readRDS(obj_save_to)
    pml_da_results <- fread(da_result_save_to)
    pml_da_deg_results <- fread(da_deg_result_save_to)
  } else {
    cat("[I]: Building Milo object and estimating DA ...\n")

    # Create millo object and build up the neighborhoods
    set.seed(RANDOM_SEED)
    pml_milo <- Milo(pbmc_int_sce[, tar_cells])
    pml_milo <- buildGraph(pml_milo, k = nhood_k, d = reductdim_d, reduced.dim = pca_source)
    pml_milo <- makeNhoods(pml_milo, prop = nhood_cell_prop, k = nhood_k, d = reductdim_d, reduced_dims = pca_source, refined = TRUE, refinement_scheme="graph")
    pml_milo <- countCells(pml_milo, samples = "Sample_by_Timepoint", meta.data = as.data.frame(colData(pml_milo)))
    pml_milo <- calcNhoodDistance(pml_milo, d = reductdim_d, reduced.dim = pca_source)
    pml_milo <- buildNhoodGraph(pml_milo)
    pml_milo <- calcNhoodExpression(pml_milo)

    # Estimate differential abundance based on the design matrix
    pml_design <- colData(pml_milo) %>%
      as.data.table() %>%
      dplyr::select(Sample_by_Timepoint, Response, Treatment) %>%
      dplyr::distinct() %>%
      tibble::column_to_rownames("Sample_by_Timepoint") %>%
      (function(m) {
        if (per_comp == "Rs_vs_NRs") {
          dplyr::filter(m, Treatment == "BL") %>% dplyr::select(-Treatment)
        } else if (per_comp == "Responder_by_time") {
          dplyr::filter(m, Response == "Responder") %>% dplyr::select(-Response)
        }
      })
    pml_fml <- colnames(pml_design) %>% paste(collapse = " + ") %>% paste("~", .) %>% formula()
    pml_da_results <- testNhoods(pml_milo, design = pml_fml, design.df = pml_design, reduced.dim = pca_source)
    pml_da_results <- tryCatch(annotateNhoods(pml_milo, pml_da_results, coldata_col = "predicted.celltype.l1"), error = function(e) {cat(e$message, "\n"); return(pml_da_results)})
    pml_da_results <- tryCatch(groupNhoods(pml_milo, pml_da_results, da.fdr = alpha), error = function(e) {cat(e$message, "\n"); return(pml_da_results)})

    # Test the differential expression
    tar_nhoods <- pml_da_results %>% dplyr::mutate(tar_nhoods = FDR < alpha) %>% dplyr::pull(tar_nhoods)
    tar_features <- rownames(pml_milo) %>% purrr::keep(~ ! .x %in% chrom_xy_genes)
    pml_da_deg_results <- testDiffExp(pml_milo, pml_da_results, design = pml_fml, subset.row = tar_features, meta.data = data.frame(colData(pml_milo)), subset.nhoods = tar_nhoods) %>%
      lapply(function(x) x %>% tibble::rownames_to_column("gene_symbol")) %>%
      Reduce(rbind, .) %>%
      dplyr::mutate(comparison = per_comp)

    # Save the Milo object
    saveRDS(pml_milo, obj_save_to)

    # Save the differential abundance results
    fwrite(pml_da_results, da_result_save_to)

    # save the differentially expressed results of DA neighborhoods
    fwrite(pml_da_deg_results, da_deg_result_save_to)
  }

  # Marker features from DA
  # KLRG1+, IL7R-, short-lived effector CD8+ T cells
  # CD8+ TEM, lack CCR7 and CD62L (CD62L is encoded by SELL)
  mark_features <- pml_da_deg_results %>% dplyr::filter(abs(logFC) > 0.45, adj.P.Val < 0.05, AveExpr > 1) %>%
    dplyr::pull(gene_symbol) %>% c(important_markers, "SYNE1", "CST7", "HCST", "ITGB1", "IL32", "TRBC2", "CD3G", "CD2") %>%
    unique %>% purrr::keep(~ .x %in% rownames(pml_milo))
  highlight_features <- c("SYNE1", "CST7", "HCST", "ITGB1", "IL32", "TRBC2", "CD3G", "CD2", important_markers)

  p <- plotNhoodExpressionDA(pml_milo, pml_da_results, features = mark_features, show_rownames = FALSE, scale_to_1 = TRUE, cluster_features = TRUE, highlight_features = highlight_features) +
    scale_fill_gradient(name = "Ave Expr", low = rbgc[5], high = rbgc[11]) &
    theme(legend.position = "right", legend.box.margin = margin(0, 0, 0, 40), legend.box.just = "right")
  save_to <- file.path(plot_dir, paste0(save_token, ".nh_marker_gene_expression.v3.pdf"))
  ggsave(save_to, plot = p, width = 10, height = 14)

  # Show distribution of neighborhood sizes
  save_to <- file.path(plot_dir, paste0(save_token, ".nhood_size_hist.pdf"))
  p <- plotNhoodSizeHist(pml_milo)
  ggsave(save_to, plot = p, width = 5, height = 5)

  # Show distribution of p-values of abundance analysis for each neighborhood
  save_to <- file.path(plot_dir, paste0(save_token, ".p_value_distribution.pdf"))
  p <- ggplot(pml_da_results, aes(PValue)) + geom_histogram(bins=50)
  ggsave(save_to, plot = p, width = 5, height = 5)

  # Show neighborhood graph by response
  save_to <- file.path(plot_dir, paste0(save_token, ".nh_graph.by_response.v2.pdf"))
  umap_pl <- plotReducedDim(pml_milo, dimred = umap_source, colour_by = umap_color_by, text_by = "predicted.celltype.l1", text_size = 5, point_size = 2) + scale_color_manual(values = c("Responder" = npg3[1], "Nonresponder" = npg3[2])) + guides(fill="none")
  nh_graph_pl <- plotNhoodGraphDA(pml_milo, pml_da_results, layout = umap_source, alpha = alpha) + scale_fill_gradient2(low = rbgc[11], mid = rbgc[6], high = rbgc[1]) + labs(fill = "Log2FC")
  p <- umap_pl + nh_graph_pl + plot_layout(guides = "collect") & theme(legend.text = element_text(size = 10), legend.title = element_text(size = 12))
  ggsave(save_to, plot = p, width = 12, height = 6)

  # Show neighborhood groups by response
  if ("NhoodGroup" %in% colnames(pml_da_results)) {
    save_to <- file.path(plot_dir, paste0(save_token, ".nh_graph.group_nhoods.by_response.pdf"))
    umap_pl <- plotReducedDim(pml_milo, dimred = umap_source, colour_by = umap_color_by, text_by = "predicted.celltype.l1", text_size = 5, point_size = 1) + guides(fill="none")
    nh_graph_pl <- plotNhoodGroups(pml_milo, pml_da_results, layout = umap_source, alpha = alpha)
    p <- umap_pl + nh_graph_pl + plot_layout(guides="collect")
    ggsave(save_to, plot = p, width = 12, height = 6)
  }

  # Show neighborhood graph by predicted cell type
  save_to <- file.path(plot_dir, paste0(save_token, ".nh_graph.by_response.predicted_celltype_l2.pdf"))
  umap_pl <- plotReducedDim(pml_milo, dimred = umap_source, colour_by = "predicted.celltype.l2", text_by = "predicted.celltype.l2", text_size = 5, point_size = 1) + guides(fill="none")
  nh_graph_pl <- plotNhoodGraphDA(pml_milo, pml_da_results, layout = umap_source, alpha = alpha)
  p <- umap_pl + nh_graph_pl + plot_layout(guides="collect")
  ggsave(save_to, plot = p, width = 12, height = 6)

  # Show neighborhood graph by co-factors
  sample_match <- c("PML0022" = "PML-1", "PML0017" = "PML-2", "PML0025" = "PML-3", "PML0002" = "PML-4", "PML0009" = "PML-5")
  new_patient_id <- sample_match[pml_milo$Vireo_assignment]
  names(new_patient_id) <- NULL
  pml_milo$patient_id <- new_patient_id

  save_to <- file.path(plot_dir, paste0(save_token, ".nh_graph.by_co_factors.pdf"))
  umap_id_pl <- plotReducedDim(pml_milo, dimred = umap_source, colour_by = "patient_id", text_by = "predicted.celltype.l1", text_size = 4, point_size = 0.25) + scale_color_npg() + labs(color = "Patient", title = "Patient")
  umap_age_pl <- plotReducedDim(pml_milo, dimred = umap_source, colour_by = "Age", text_by = "predicted.celltype.l1", text_size = 4, point_size = 0.25) + labs(color = "Age", title = "Age")
  umap_sex_pl <- plotReducedDim(pml_milo, dimred = umap_source, colour_by = "Sex", text_by = "predicted.celltype.l1", text_size = 4, point_size = 0.25) + scale_color_jama() + labs(color = "Gender", title = "Gender")
  umap_pool_pl <- plotReducedDim(pml_milo, dimred = umap_source, colour_by = "SequencingPool", text_by = "predicted.celltype.l1", text_size = 4, point_size = 0.25) + scale_color_lancet() + labs(color = "Library", title = "Library")
  nh_graph_pl <- plotNhoodGraphDA(pml_milo, pml_da_results, layout = umap_source, alpha = alpha) 
  p <- (umap_id_pl + umap_pool_pl + umap_age_pl + umap_sex_pl + plot_layout(ncol = 2, guides = "collect")) / nh_graph_pl + plot_layout(ncol = 2)
  ggsave(save_to, plot = p, width = 12, height = 5.5)

  # Show differential abundance beeswarm of neightborhood by predicted cell type
  save_to <- file.path(plot_dir, paste0(save_token, ".da_beeswarm.pdf"))
  if (min(pml_da_results$SpatialFDR) <= alpha) {
    p <- plotDAbeeswarm(pml_da_results, group.by = "predicted.celltype.l1", alpha = alpha) + labs(x = NULL) + theme(axis.title = element_text(size = 20))
    ggsave(save_to, plot = p, width = 7.5, height = 5)
  }
}


#
## Check the identified DA neighborhoods
#
# per_comp <- "Responder_by_time"
per_comp <- "Rs_vs_NRs"
sel_pml_da_results <- pml_da_results %>% dplyr::filter(SpatialFDR <= alpha)
sel_nhood_tab <- sel_pml_da_results %>%
  dplyr::pull(Nhood) %>%
  data.frame(Nhood = ., Nhood_idx = unlist(pml_milo@nhoodIndex[.])) %>%
  dplyr::mutate(Nhood_idx = paste0("Nhood_idx_", Nhood_idx))

da_nhoods <- pml_milo@nhoods %>%
  as.data.frame() %>%
  dplyr::rename_with(~paste0("Nhood_idx_", .x)) %>%
  dplyr::mutate(Cellbarcode = rownames(.)) %>%
  dplyr::select(dplyr::one_of(c("Cellbarcode", sel_nhood_tab$Nhood_idx))) %>%
  tidyr::pivot_longer(dplyr::starts_with("Nhood_idx_"), names_to = "Nhood_idx", values_to = "BelongTo") %>%
  dplyr::filter(BelongTo == 1) %>%
  dplyr::left_join(sel_nhood_tab, by = "Nhood_idx") %>%
  dplyr::left_join(sel_pml_da_results, by = "Nhood") %>%
  dplyr::select(-c(BelongTo, Nhood_idx, predicted.celltype.l1, predicted.celltype.l1_fraction))

# pbmc_int@meta.data %>% as.data.frame() %>% dplyr::select(predicted.celltype.l2, Response) %>% table()
pbmc_meta_data <- pbmc_int@meta.data %>% as.data.frame() %>% dplyr::mutate(Cellbarcode = rownames(.)) %>% dplyr::inner_join(da_nhoods, by = "Cellbarcode")

tar_cells <- dplyr::left_join(da_nhoods, pbmc_meta_data, by = "Cellbarcode") %>% dplyr::filter(predicted.celltype.l1 == "CD8 T") %>% dplyr::pull(Cellbarcode) %>% unique()
p <- DimPlot(pbmc_int[, tar_cells], group.by = "predicted.celltype.l2", label = TRUE, pt.size = 0.5)
save_to <- file.path(plot_dir, paste0(save_token, ".da_nhoods.l2_response.pdf"))
ggsave(save_to, plot = p, width = 7, height = 5)

sample_match <- c("PML0022" = "PML-1", "PML0017" = "PML-2", "PML0025" = "PML-3", "PML0002" = "PML-4", "PML0009" = "PML-5")
pbmc_int$patient_id <- sample_match[pbmc_int$Vireo_assignment] %>% as.character()
p <- DimPlot(pbmc_int[, tar_cells], group.by = "patient_id", label = TRUE, pt.size = 0.5)
save_to <- file.path(plot_dir, paste0(save_token, ".da_nhoods.colored_by_donors.v2.pdf"))
ggsave(save_to, plot = p, width = 6, height = 5)

if (per_comp == "Rs_vs_NRs" && per_omics == "pml_citeseq") {
  # DE genes between Responder vs Non-responder in CD8 TEM cells
  tar_timepoint <- {if (per_comp == "Rs_vs_NRs") "BL" else c("3M", "6W")}
  tp_token <- paste0(tar_timepoint, collapse = "_")
  detab_save_to <- file.path(proj_dir, paste0("outputs/analysis/differential_abundance/pbmc.pml_citeseq.", tp_token, ".de_gene.CD8_TEM.response_vs_nonresponse.csv"))
  if (file.exists(detab_save_to)) {
    deg_tem_tab <- fread(detab_save_to)
  } else {
    tar_cells <- pbmc_int@meta.data %>% dplyr::filter(predicted.celltype.l2 == "CD8 TEM", predicted.celltype.l1 == "CD8 T", Timepoint %in% tar_timepoint) %>% rownames()
    sub_pbmc_int <- pbmc_int[, tar_cells]

    deg_tem_tab <- FindMarkers(sub_pbmc_int, ident.1 = "Responder", ident.2 = "Non-responder", group.by = "Response", recorrect_umi = FALSE, logfc.threshold = 0) %>%
      dplyr::mutate(gene_symbol = rownames(.), p_val_adj = p.adjust(p_val, method = "bonferroni")) %>%
      dplyr::filter(!gene_symbol %in% chrom_xy_genes) %>%
      as.data.table()
    fwrite(deg_tem_tab, detab_save_to)
  }

  # DEGs
  deg_tem_plot_tab <- dplyr::filter(deg_tem_tab, pct.1 > 0.01, pct.2 > 0.01)
  vp_plot_save_to <- file.path(plot_dir, paste0("daa.pml_citeseq.", tp_token, ".de_gene.CD8_TEM.response_vs_nonresponse.vocano_plot.v3.pdf"))
  p_vocano <- deg_tem_plot_tab %>% ggplot() +
    geom_point(mapping = aes(x = avg_log2FC, y = -log10(p_val_adj)), data = deg_tem_plot_tab, size = 0.25) +
    geom_text_repel(mapping = aes(x = avg_log2FC, y = -log10(p_val_adj), label = gene_symbol), data = dplyr::slice_min(deg_tem_plot_tab, p_val_adj, n = 15), box.padding = 0.5, min.segment.length = 0, max.overlaps = Inf)
  #if (per_comp == "Rs_vs_NRs") {
  #    p_vocano <- p_vocano +
  #      geom_point(mapping = aes(x = avg_log2FC, y = -log10(p_val_adj)), data = dplyr::filter(deg_tem_plot_tab, gene_symbol == "IFNG"), size = 2, color = "red") +
  #      geom_text_repel(mapping = aes(x = avg_log2FC, y = -log10(p_val_adj), label = gene_symbol), data = dplyr::filter(deg_tem_plot_tab, gene_symbol == "IFNG"), color = "red", min.segment.length = 0)
  #}
  p_vocano <- p_vocano + geom_vline(xintercept = c(0.25, -0.25), linetype = "dashed") + geom_hline(yintercept = -log10(0.05), linetype = "dashed") + theme_classic()
  ggsave(vp_plot_save_to, plot = p_vocano, width = 4, height = 5)

  # Differnetial expression between responder and non-responder of IFNG in CD8 TEM
  vd_plot_save_to <- file.path(plot_dir, "daa.pml_citeseq.BL.CD8_TEM.ifng_genes.vln_plot.pdf")
  p_vln <- VlnPlot(sub_pbmc_int, assay = "RNA", features = "IFNG", group.by = "Response", pt.size = 0, adjust = 1.5) +
    ylim(c(-0.01, 4.5)) +
    geom_signif(annotation = formatC(dplyr::filter(deg_tem_plot_tab, gene_symbol == "IFNG")$p_val_adj, digits = 1), y_position = 4.4, xmin = 1, xmax = 2) +
    theme(legend.title = element_text(size = 12), axis.text.x = element_blank(), axis.ticks.x = element_blank(), plot.margin = margin(0, 0, 0, 0)) +
    labs(fill = "Response", title = NULL, x = NULL) +
    scale_fill_jama()

  p_dot <- DotPlot(sub_pbmc_int, features = "IFNG", group.by = "Response", cols = "Spectral") +
    labs(y = NULL, x = NULL) +
    theme(plot.margin = margin(0, 0, 0, 0), legend.title = element_text(size = 12), axis.line.y.left = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank()) +
    scale_size(name = "Perc. exp.", range = c(3, 9)) +
    scale_y_discrete(breaks = c("Responder", "Non-responder"), labels = c("Rs", "NRs")) +
    coord_flip()
  p <- p_vln + p_dot + plot_layout(guides = "collect", heights = c(7, 1))
  ggsave(vd_plot_save_to, plot = p, width = 4, height = 5)

  # Functional enrichment analysis of DEGs
  sig_deg_tem_tab <- deg_tem_tab %>% dplyr::filter(p_val_adj < 0.05, pct.1 > 0.05, pct.2 > 0.05) %>% dplyr::mutate(Direction = dplyr::if_else(avg_log2FC < 0, "dw", "up"))
  sig_deg_tem_tab$Direction %>% table()
  for (pp in c("up", "dw")) {
    bpgo_tar_genes <- dplyr::filter(sig_deg_tem_tab, Direction == pp) %>% dplyr::pull(gene_symbol)
    bp_ego <- enrichGO(gene = bpgo_tar_genes, ont = "BP", keyType = "SYMBOL", OrgDb = "org.Hs.eg.db", readable = TRUE, pvalueCutoff = 0.05, qvalueCutoff = 0.05)

    if (dim(bp_ego)[1] != 0) {
      gp_save_to <- file.path(plot_dir, paste("daa.pml_citeseq", tp_token, "de_gene.CD8_TEM.response_vs_nonresponse", pp, "go_enrichment.csv", sep = "."))
      bp_ego@result %>% fwrite(gp_save_to)

      bp_ego_2 <- pairwise_termsim(bp_ego)
      tp_p <- treeplot(bp_ego_2, showCategory = 25, offset = 20, offset_tiplab = 1) +
        scale_color_gradient(low = rbgc[1], high = "gray90") +
        labs(color = "Adj. P value", sizes = "Nr. of genes") +
        theme(legend.position = "left")
      tp_save_to <- file.path(plot_dir, paste("daa.pml_citeseq", tp_token, "de_gene.CD8_TEM.response_vs_nonresponse", pp, "go_enrichment.treeplot.pdf", sep = "."))
      ggsave(tp_save_to, plot = tp_p, height = 6, width = 13.5)

      go_tar_genes <- dplyr::filter(sig_deg_tem_tab, Direction == pp) %>% dplyr::arrange(desc(avg_log2FC)) %>% dplyr::pull(avg_log2FC, gene_symbol)
      hp_save_to <- file.path(plot_dir, paste("daa.pml_citeseq", tp_token, "de_gene.CD8_TEM.response_vs_nonresponse", pp, "go_enrichment.heatmap.pdf", sep = "."))
      hp_p <- heatplot(bp_ego, showCategory = 25, foldChange = go_tar_genes)
      ggsave(hp_save_to, plot = hp_p, height = 5, width = 12.5)
    }

    bp_ego_save_to <- file.path(object_dir, paste("daa.pml_citeseq", tp_token, "de_gene.CD8_TEM.response_vs_nonresponse", pp, "go_enrichment.rds", sep = "."))
    saveRDS(bp_ego, bp_ego_save_to)
  }

  # Gene set enrichment, GO
  for (pgo in c("BP", "CC", "MF")) {
    bp_gsa_go_save_to <- file.path(plot_dir, "GSEA", paste("daa.pml_citeseq", tp_token, "de_gene.CD8_TEM.response_vs_nonresponse.go_gsea", pgo, "rds", sep = "."))
    if (!file.exists(bp_gsa_go_save_to)) {
      gse_tar_genes <- sig_deg_tem_tab %>% dplyr::pull(avg_log2FC, gene_symbol) %>% sort(decreasing = TRUE)
      bp_gsa_go <- gseGO(gse_tar_genes, OrgDb = "org.Hs.eg.db", ont = pgo, keyType = "SYMBOL", minGSSize = 15, maxGSSize = 500, pvalueCutoff = 0.05)
      saveRDS(bp_gsa_go, bp_gsa_go_save_to)
    } else {
      bp_gsa_go <- readRDS(bp_gsa_go_save_to)
    }

    # GSEA plot of NES
    nes_tab <- bp_gsa_go@result %>% as.data.table() %>%
      dplyr::mutate(n_genes = .$core_enrichment %>% purrr::map(~stringr::str_split(.x, "/", simplify=T) %>% length) %>% unlist()) %>%
      dplyr::mutate(GeneRatio = n_genes / setSize, GeneRatioRank = rank(GeneRatio)) %>%
      dplyr::mutate(Description = forcats::fct_reorder(Description, -log10(p.adjust)))

    selected_terms <- nes_tab %>% dplyr::filter(NES > 0)
    p <- ggplot(nes_tab) +
      geom_bar(aes(x = -log10(p.adjust), y = Description, fill = NES, color=NES), stat = "identity") +
      geom_text_repel(aes(x = -log10(p.adjust), y = Description, label = ID), nudge_x = 1.5, hjust = 0, direction = "y", min.segment.length = 0, data = selected_terms) +
      scale_fill_gradient2(low = rbgc[11], mid=rbgc[6], high=rbgc[1]) +
      scale_color_gradient2(low = rbgc[11], mid=rbgc[6], high=rbgc[1]) +
      labs(x = "-log10(Adjusted-p-value)", y = paste0("Gene ontology term (", pgo, ")")) +
      theme_classic() +
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
    p_save_to <- file.path(plot_dir, "GSEA", paste("daa.pml_citeseq", tp_token, "de_gene.CD8_TEM.response_vs_nonresponse.go_gsea", pgo, "NES.pdf", sep = "."))
    ggsave(p_save_to, plot = p, height = 6, width = 4.5)

    # # GSEA plot per GO term
    # for (per_idx in 1:length(bp_gsa_go$Description)) {
    #   per_desc <- bp_gsa_go$Description[per_idx]
    #   per_desc_name <- stringr::str_replace_all(per_desc, " ", "_")
    #   gsa_p_save_to <- file.path(plot_dir, "GSEA", paste("daa.pml_citeseq", tp_token, "de_gene.CD8_TEM.response_vs_nonresponse.go_gsea", pgo, per_desc_name, "pdf", sep = "."))
    #   gsa_p <- gseaplot2(bp_gsa_go, geneSetID = per_idx, title = per_desc, pvalue_table = TRUE)
    #   ggsave(gsa_p_save_to, plot = gsa_p, height = 4.5, width = 8)
    # }
  }

  # Selecte gene set enrichment, BP
  if (per_comp == "Responder_by_time") {
    cat("TODO:\n")
  } else if (per_comp == "Rs_vs_NRs") {
    bp_gsa_go_save_to <- file.path(plot_dir, "GSEA", paste0("daa.pml_citeseq.", tp_token, ".de_gene.CD8_TEM.response_vs_nonresponse.go_gsea.BP.rds"))
    bp_gsa_go <- readRDS(bp_gsa_go_save_to)
    selected_gp_desc <- c("positive regulation of response to stimulus", "positive regulation of apoptotic process", "positive regulation of cell migration")
    selected_gp_idx <- which(bp_gsa_go$Description %in% selected_gp_desc)
    selected_gp_desc_short <- bp_gsa_go$Description[selected_gp_idx] %>% stringr::str_replace_all("positive regulation", "Pos. reg.") %>% purrr::set_names(bp_gsa_go$Description[selected_gp_idx])
    gsa_p_save_to <- file.path(plot_dir, "GSEA", paste0("daa.pml_citeseq.BL.de_gene.CD8_TEM.response_vs_nonresponse.go_gsea.pdf"))
    gsa_p <- gseaplot2(bp_gsa_go, geneSetID = selected_gp_idx) & scale_color_manual(values = lacet3, labels = selected_gp_desc_short) & theme(legend.text = element_text(size = 11))
    ggsave(gsa_p_save_to, plot = gsa_p, height = 4.5, width = 6.75)

    gsa_p_save_to <- file.path(plot_dir, "GSEA", paste0("daa.pml_citeseq.BL.de_gene.CD8_TEM.response_vs_nonresponse.go_gsea.dot.pdf"))
    p <- dotplot(bp_gsa_go, showCategory = 20)
    ggsave(gsa_p_save_to, plot = p, height = 8.5, width = 6.75)
 
    # GSEA, DE genes at baseline between responders and non-responders
    deg_tem_tab <- file.path(proj_dir, "outputs/analysis/differential_abundance", paste0("pbmc.pml_citeseq.", tp_token, ".de_gene.CD8_TEM.response_vs_nonresponse.csv")) %>%
      fread() %>% dplyr::filter(!gene_symbol %in% chrom_xy_genes, p_val_adj < 0.05, abs(avg_log2FC) > 0.1)
    deg_all_tab <- file.path(proj_dir, "outputs/analysis/deg/tables/pbmc.pml_citeseq.integrated.de_gene.csv") %>%
      fread() %>% dplyr::filter(donors == "All", !gene_symbol %in% chrom_xy_genes, comparison == "Rs.vs.NRs_BL", p_val_adj < 0.05, abs(avg_log2FC) > 0.1)

    reference_gene_sets <- deg_all_tab %>%
      dplyr::group_by(celltype) %>%
      # dplyr::group_by(celltype, direction) %>%
      dplyr::summarize(gene_set = list(gene_symbol)) %>%
      dplyr::mutate(gene_set_name = stringr::str_remove_all(celltype, " ")) %>%
      # dplyr::mutate(gene_set_name = paste0(stringr::str_remove_all(celltype, " "), "_", direction)) %>%
      dplyr::pull(gene_set, gene_set_name)

    deg_tem <- deg_tem_tab %>% dplyr::pull(avg_log2FC, gene_symbol) %>% sort(decreasing = TRUE)

    p <- NULL
    for (pct in c("CD8T", "CD4T", "Mono", "NK", "B")) {
      if (is.null(p)) {
        p <- plotEnrichment(reference_gene_sets[[pct]], deg_tem) + labs(x = NULL, y = paste0(pct, " (ES)"))
      } else {
        p <- p / plotEnrichment(reference_gene_sets[[pct]], deg_tem) + labs(x = NULL, y = paste0(pct, " (ES)"))
      }
    }

    p <- p + labs(x = "Ranks")
    tem_gsea_save_to <- file.path(plot_dir, "GSEA", paste0("daa.pml_citeseq.CD8_TEM.response_vs_nonresponse.gsea_in_BL_deg.pdf"))
    ggsave(tem_gsea_save_to, plot = p, height = 8, width = 6)
  }
}
