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
  library(clusterProfiler)

  library(miloR)
  library(SingleCellExperiment)
  library(Seurat)
  library(scran)
  library(scater)
})


rm(list = ls()); gc()
overwrite <- FALSE
proj_dir <- "~/Documents/projects/wp_pml"
plot_dir <- file.path(proj_dir, "outputs/analysis/differential_abundance/plots")
object_dir <- file.path(proj_dir, "outputs/analysis/differential_abundance/objects")
set.seed(3141592)

my_colors <- colorRampPalette(c("red"))(10)


# per_omics <- "pml_rnaseq"; pca_source <- "PCA"; umap_source <- "UMAP"
# per_omics <- "pml_rnaseq"; pca_source <- "REF.SPCA"; umap_source <- "REF.UMAP"
# per_omics <- "pml_citeseq"; pca_source <- "REF.SPCA"; umap_source <- "REF.UMAP"
per_omics <- "pml_citeseq"; pca_source <- "INT_PCA"; umap_source <- "WNN.UMAP"
pca_source_dict <- c("PCA" = "internal_pca", "INT_PCA" = "internal_pca", "REF.SPCA" = "external_pca")

pbmc_int <- file.path(proj_dir, "outputs/analysis/integrated/objects", paste0(per_omics, ".rds")) %>% readRDS()
pbmc_int_sce <- as.SingleCellExperiment(pbmc_int, assay = "RNA")
# pbmc_int_sce <- as.SingleCellExperiment(pbmc_int, assay = "SCT")
pbmc_int_sce$Response <- factor(pbmc_int_sce$Response, levels = c("Non-responder", "Responder"))

save_to <- file.path(plot_dir, paste("daa", per_omics, "ref_umap.scran.pdf", sep = "."))
p <- plotReducedDim(pbmc_int_sce, "REF.UMAP", colour_by = "Response")
ggsave(save_to, plot = p, width = 5, height = 5)

save_to <- file.path(plot_dir, paste("daa", per_omics, "umap.scran.pdf", sep = "."))
p <- plotReducedDim(pbmc_int_sce, umap_source, colour_by = "Response")
ggsave(save_to, plot = p, width = 5, height = 5)


# Per time points by response
alpha <- 0.1
nhood_k <- 50
reductdim_d <- 30
for (ptp in c("BL")) {
  save_token <- paste("daa", per_omics, ptp, pca_source_dict[pca_source], sep = ".")
  tar_cells <- pbmc_int_sce %>% colData() %>% as.data.frame() %>% dplyr::filter(Timepoint == ptp) %>% rownames()

  obj_save_to <- file.path(object_dir, paste0(save_token, ".differential_abundance_test_object.rds"))
  da_result_save_to <- file.path(object_dir, paste0(save_token, ".differential_abundance_test_results.csv"))
  if (!file.exists(obj_save_to) || !file.exists(da_result_save_to) || overwrite) {
    cat("[I]: Creating from the raw data ...\n")
    pml_milo <- Milo(pbmc_int_sce[, tar_cells])
    pml_milo <- buildGraph(pml_milo, k = nhood_k, d = reductdim_d, reduced.dim = pca_source)
    pml_milo <- makeNhoods(pml_milo, prop = 0.5, k = nhood_k, d = reductdim_d, refined = TRUE, reduced_dims = pca_source)
    pml_milo <- countCells(pml_milo, meta.data = as.data.frame(colData(pml_milo)), sample = "Vireo_assignment")
    pml_milo <- calcNhoodDistance(pml_milo, d = reductdim_d, reduced.dim = pca_source)
    pml_milo <- buildNhoodGraph(pml_milo)
    saveRDS(pml_milo, obj_save_to)

    pml_design <- colData(pml_milo) %>%
      as.data.table() %>%
      dplyr::filter(Timepoint == ptp) %>%
      dplyr::select(Vireo_assignment, Response) %>%
      dplyr::distinct() %>%
      as.data.frame() %>%
      dplyr::mutate(Response = factor(Response, levels = c("Non-responder", "Responder"))) %>%
      tibble::column_to_rownames("Vireo_assignment")
    pml_da_results <- testNhoods(pml_milo, design = ~ Response, design.df = pml_design, reduced.dim = pca_source)
    pml_da_results <- annotateNhoods(pml_milo, pml_da_results, coldata_col = "predicted.celltype.l1")
    tryCatch(pml_da_results <- groupNhoods(pml_milo, pml_da_results, max.lfc.delta = 1), error = function(e) cat(e$message, "\n"))
    fwrite(pml_da_results, da_result_save_to)

    save_to <- file.path(plot_dir, paste0(save_token, ".p_value_distribution.pdf"))
    p <- ggplot(pml_da_results, aes(PValue)) + geom_histogram(bins=50)
    ggsave(save_to, plot = p, width = 5, height = 5)

    save_to <- file.path(plot_dir, paste0(save_token, ".nhood_size_hist.pdf"))
    p <- plotNhoodSizeHist(pml_milo)
    ggsave(save_to, plot = p, width = 5, height = 5)

    if ("NhoodGroup" %in% colnames(pml_da_results)) {
      save_to <- file.path(plot_dir, paste0(save_token, ".nh_graph.group_nhoods.by_response.pdf"))
      umap_pl <- plotReducedDim(pml_milo, dimred = umap_source, colour_by = "Response", text_by = "predicted.celltype.l1", text_size = 5, point_size = 1) + guides(fill="none")
      nh_graph_pl <- plotNhoodGroups(pml_milo, pml_da_results, layout = umap_source, alpha = alpha)
      p <- umap_pl + nh_graph_pl + plot_layout(guides="collect")
      ggsave(save_to, plot = p, width = 12, height = 6)
    }

    save_to <- file.path(plot_dir, paste0(save_token, ".nh_graph.by_response.pdf"))
    umap_pl <- plotReducedDim(pml_milo, dimred = umap_source, colour_by = "Response", text_by = "predicted.celltype.l1", text_size = 5, point_size = 1) + guides(fill="none")
    nh_graph_pl <- plotNhoodGraphDA(pml_milo, pml_da_results, layout = umap_source, alpha = alpha)
    p <- umap_pl + nh_graph_pl + plot_layout(guides="collect")
    ggsave(save_to, plot = p, width = 12, height = 6)

    save_to <- file.path(plot_dir, paste0(save_token, ".nh_graph.by_response.predicted_celltype_l2.pdf"))
    umap_pl <- plotReducedDim(pml_milo, dimred = umap_source, colour_by = "predicted.celltype.l2", text_by = "predicted.celltype.l2", text_size = 5, point_size = 1) + guides(fill="none")
    nh_graph_pl <- plotNhoodGraphDA(pml_milo, pml_da_results, layout = umap_source, alpha = alpha)
    p <- umap_pl + nh_graph_pl + plot_layout(guides="collect")
    ggsave(save_to, plot = p, width = 12, height = 6)

    save_to <- file.path(plot_dir, paste0(save_token, ".nh_graph.by_co_factors.pdf"))
    umap_id_pl <- plotReducedDim(pml_milo, dimred = umap_source, colour_by = "Vireo_assignment", text_by = "predicted.celltype.l1", text_size = 4, point_size = 0.25) + scale_color_npg() + labs(color = "Patient", title = "Patient")
    umap_age_pl <- plotReducedDim(pml_milo, dimred = umap_source, colour_by = "Age", text_by = "predicted.celltype.l1", text_size = 4, point_size = 0.25) + labs(color = "Age", title = "Age")
    umap_sex_pl <- plotReducedDim(pml_milo, dimred = umap_source, colour_by = "Sex", text_by = "predicted.celltype.l1", text_size = 4, point_size = 0.25) + scale_color_jama() + labs(color = "Gender", title = "Gender")
    umap_pool_pl <- plotReducedDim(pml_milo, dimred = umap_source, colour_by = "pool", text_by = "predicted.celltype.l1", text_size = 4, point_size = 0.25) + scale_color_lancet() + labs(color = "Library", title = "Library")
    nh_graph_pl <- plotNhoodGraphDA(pml_milo, pml_da_results, layout = umap_source, alpha = alpha) 
    p <- (umap_id_pl + umap_pool_pl + umap_age_pl + umap_sex_pl + plot_layout(ncol = 2, guides = "collect")) / nh_graph_pl + plot_layout(ncol = 2)
    ggsave(save_to, plot = p, width = 12, height = 5.5)

    save_to <- file.path(plot_dir, paste0(save_token, ".da_beeswarm.pdf"))
    if (min(pml_da_results$SpatialFDR) <= alpha) {
      p <- plotDAbeeswarm(pml_da_results, group.by = "predicted.celltype.l1", alpha = alpha)
      ggsave(save_to, plot = p, width = 6, height = 5)
    }
  } else {
    cat("[I]: Loading from the disk ...\n")
    pml_milo <- readRDS(obj_save_to)
    pml_da_results <- fread(da_result_save_to) %>% dplyr::mutate(NhoodGroup = as.factor(NhoodGroup))
  }
}

# Check the identified CD8 TEM
sel_pml_da_results <- pml_da_results %>% dplyr::filter(SpatialFDR <= 0.05)
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
pbmc_meta_data <- pbmc_int@meta.data %>%
  as.data.frame() %>%
  dplyr::mutate(Cellbarcode = rownames(.)) %>%
  dplyr::inner_join(da_nhoods, by = "Cellbarcode")

tar_cells <- da_nhoods %>% dplyr::pull(Cellbarcode) %>% unique()
p <- DimPlot(pbmc_int[, tar_cells], group.by = "predicted.celltype.l2", label = TRUE)
save_to <- file.path(plot_dir, paste0(save_token, ".da_nhoods.l2_response.pdf"))
ggsave(save_to, plot = p, width = 6, height = 5)

p <- DimPlot(pbmc_int[, tar_cells], group.by = "Vireo_assignment", label = TRUE)
save_to <- file.path(plot_dir, paste0(save_token, ".da_nhoods.colored_by_donors.pdf"))
ggsave(save_to, plot = p, width = 6, height = 5)


# DE genes between Responder vs Non-responder in CD8 TEM cells
detab_save_to <- file.path(proj_dir, "outputs/analysis/differential_abundance/pbmc.pml_citeseq.BL.de_gene.CD8_TEM.response_vs_nonresponse.csv")
if (file.exists(detab_save_to)) {
  deg_tem_tab <- fread(detab_save_to)
} else {
  tar_cells <- pbmc_int@meta.data %>% dplyr::filter(predicted.celltype.l2 == "CD8 TEM", predicted.celltype.l1 == "CD8 T", Timepoint == "BL") %>% rownames()
  sub_pbmc_int <- pbmc_int[, tar_cells]

  deg_tem_tab <- FindMarkers(sub_pbmc_int, ident.1 = "Responder", ident.2 = "Non-responder", group.by = "Response", recorrect_umi = FALSE, test.use = "t", min.pct = 0.10, logfc.threshold = 0) %>%
    dplyr::mutate(gene_symbol = rownames(.), p_val_adj = p.adjust(p_val, method = "bonferroni")) %>%
    as.data.table()
  fwrite(deg_tem_tab, detab_save_to)
}

# DEGs
vp_plot_save_to <- file.path(plot_dir, "daa.pml_citeseq.BL.de_gene.CD8_TEM.response_vs_nonresponse.vocano_plot.pdf")
p_vocano <- ggplot() +
  geom_point(mapping = aes(x = avg_log2FC, y = -log10(p_val_adj)), data = deg_tem_tab, size = 0.25) +
  geom_text_repel(mapping = aes(x = avg_log2FC, y = -log10(p_val_adj), label = gene_symbol), data = dplyr::slice_min(deg_tem_tab, p_val_adj, n = 15), box.padding = 0.5, min.segment.length = 0, max.overlaps = Inf) +
  geom_point(mapping = aes(x = avg_log2FC, y = -log10(p_val_adj)), data = dplyr::filter(deg_tem_tab, gene_symbol == "IFNG"), size = 2, color = "red") +
  geom_text_repel(mapping = aes(x = avg_log2FC, y = -log10(p_val_adj), label = gene_symbol), data = dplyr::filter(deg_tem_tab, gene_symbol == "IFNG"), color = "red", min.segment.length = 0) +
  geom_vline(xintercept = c(0.25, -0.25), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_classic()
ggsave(vp_plot_save_to, plot = p_vocano, width = 4, height = 5)


# Differnetial expression between responder and non-responder of IFNG in CD8 TEM
vd_plot_save_to <- file.path(plot_dir, "daa.pml_citeseq.BL.CD8_TEM.ifng_genes.vln_plot.pdf")
p_vln <- VlnPlot(sub_pbmc_int, assay = "RNA", features = "IFNG", group.by = "Response", pt.size = 0.5) +
  ylim(c(-0.01, 4.5)) +
  geom_signif(annotation = formatC(dplyr::filter(deg_tem_tab, gene_symbol == "IFNG")$p_val_adj, digits = 1), y_position = 4.4, xmin = 1, xmax = 2) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), plot.margin = margin(0, 0, 0, 0)) +
  labs(fill = "Response", title = NULL, x = NULL) +
  scale_fill_jama()
p_dot <- DotPlot(sub_pbmc_int, features = "IFNG", group.by = "Response", cols = "Spectral") +
  labs(y = NULL, x = NULL) +
  theme(plot.margin = margin(0, 0, 0, 0), axis.line.y.left = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank()) +
  scale_size(range = c(2, 6)) +
  coord_flip()
p <- p_vln + p_dot + plot_layout(guides = "collect", heights = c(7, 1))
ggsave(vd_plot_save_to, plot = p, width = 5.5, height = 4)


# Functional enrichment analysis of DEGs
sig_deg_tem_tab <- deg_tem_tab %>%
  dplyr::filter(p_val_adj < 0.05, abs(avg_log2FC) > 0.2) %>%
  dplyr::mutate(Direction = dplyr::if_else(avg_log2FC < 0, "dw", "up"))
sig_deg_tem_tab$Direction %>% table()

for (pp in c("up", "dw")) {
  bpgo_tar_genes <- dplyr::filter(sig_deg_tem_tab, Direction == pp) %>% dplyr::pull(gene_symbol)
  bp_ego <- enrichGO(gene = bpgo_tar_genes, ont = "BP", keyType = "SYMBOL", OrgDb = "org.Hs.eg.db", readable = TRUE, pvalueCutoff = 0.05, qvalueCutoff = 0.05)

  if (dim(bp_ego)[1] != 0) {
    gp_save_to <- file.path(plot_dir, paste("daa.pml_citeseq.BL.de_gene.CD8_TEM.response_vs_nonresponse", pp, "go_enrichment.csv", sep = "."))
    bp_ego@result %>% fwrite(gp_save_to)

    bp_ego_2 <- pairwise_termsim(bp_ego)
    tp_p <- treeplot(bp_ego_2, showCategory = 25, offset = 30, offset_tiplab = 1) +
      scale_color_gradient(low = "gray90", high = "red") +
      labs(color = "Adj. P value", sizes = "Nr. of genes") +
      theme(legend.position = "left")

    tp_save_to <- file.path(plot_dir, paste("daa.pml_citeseq.BL.de_gene.CD8_TEM.response_vs_nonresponse", pp, "go_enrichment.treeplot.pdf", sep = "."))
    ggsave(tp_save_to, plot = tp_p, height = 4, width = 8)

    go_tar_genes <- dplyr::filter(sig_deg_tem_tab, Direction == pp) %>% dplyr::arrange(desc(avg_log2FC)) %>% dplyr::pull(avg_log2FC, gene_symbol)
    hp_save_to <- file.path(plot_dir, paste("daa.pml_citeseq.BL.de_gene.CD8_TEM.response_vs_nonresponse", pp, "go_enrichment.heatmap.pdf", sep = "."))
    hp_p <- heatplot(bp_ego, showCategory = 20, foldChange = go_tar_genes)
    ggsave(hp_save_to, plot = hp_p, height = 3, width = 10)
  }
}


gse_tar_genes <- dplyr::arrange(sig_deg_tem_tab, desc(avg_log2FC)) %>% dplyr::pull(avg_log2FC, gene_symbol)
bp_gsa_go <- gseGO(gse_tar_genes, OrgDb = "org.Hs.eg.db", ont = "BP", keyType = "SYMBOL", minGSSize = 15, maxGSSize = 500, pvalueCutoff = 0.05)
for (per_idx in 1:length(bp_gsa_go$Description)) {
  per_desc <- bp_gsa_go$Description[per_idx]
  per_desc_name <- stringr::str_replace_all(per_desc, " ", "_")
  gsa_p_save_to <- file.path(plot_dir, "GSEA", paste0("daa.pml_citeseq.BL.de_gene.CD8_TEM.response_vs_nonresponse.go_gsea.", per_desc_name, ".pdf"))
  gsa_p <- gseaplot2(bp_gsa_go, geneSetID = per_idx, title = per_desc)
  ggsave(gsa_p_save_to, plot = gsa_p, height = 4.5, width = 7)
}

# Gene set enrichment
gsa_p_save_to <- file.path(plot_dir, "GSEA", paste0("daa.pml_citeseq.BL.de_gene.CD8_TEM.response_vs_nonresponse.go_gsea.pdf"))
gsa_p <- gseaplot2(bp_gsa_go, geneSetID = c(18, 20, 44))
ggsave(gsa_p_save_to, plot = gsa_p, height = 4.5, width = 8)

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
# [1] stats4    stats     graphics  grDevices utils     datasets  methods
# [8] base                                          
#                                               
# other attached packages:                            
#  [1] scater_1.24.0               scran_1.24.1   
#  [3] scuttle_1.6.3               sp_1.5-0   
#  [5] SeuratObject_4.1.0          Seurat_4.1.1      
#  [7] SingleCellExperiment_1.18.1 SummarizedExperiment_1.26.1
#  [9] Biobase_2.56.0              GenomicRanges_1.48.0
# [11] GenomeInfoDb_1.35.15        IRanges_2.30.0    
# [13] S4Vectors_0.34.0            BiocGenerics_0.42.0
# [15] MatrixGenerics_1.8.1        matrixStats_0.62.0      
# [17] miloR_1.4.0                 edgeR_3.38.4      
# [19] limma_3.52.4                ggsignif_0.6.4
# [21] ggsci_2.9                   dichromat_2.0-0.1  
# [23] patchwork_1.1.1             forcats_0.5.1
# [25] stringr_1.4.1               dplyr_1.0.9   
# [27] purrr_0.3.4                 readr_2.1.2      
# [29] tidyr_1.2.0                 tibble_3.1.8
# [31] ggplot2_3.4.0               tidyverse_1.3.2       
# [33] data.table_1.14.2           magrittr_2.0.3
#                                             
# loaded via a namespace (and not attached):   
#   [1] utf8_1.2.2                reticulate_1.25                                                                                                         [0/2000]
#   [3] tidyselect_1.2.0          htmlwidgets_1.5.4
#   [5] grid_4.2.0                BiocParallel_1.30.3
#   [7] Rtsne_0.16                munsell_0.5.0
#   [9] ScaledMatrix_1.4.1        codetools_0.2-18
#  [11] ica_1.0-3                 statmod_1.5.0
#  [13] future_1.26.1             miniUI_0.1.1.1
#  [15] withr_2.5.0               spatstat.random_3.0-1
#  [17] colorspace_2.0-3          progressr_0.10.1
#  [19] ROCR_1.0-11               tensor_1.5
#  [21] listenv_0.8.0             GenomeInfoDbData_1.2.8
#  [23] polyclip_1.10-0           farver_2.1.1
#  [25] parallelly_1.32.0         vctrs_0.5.1
#  [27] generics_0.1.3            R6_2.5.1
#  [29] ggbeeswarm_0.7.2          graphlayouts_0.8.0
#  [31] rsvd_1.0.5                locfit_1.5-9.6
#  [33] bitops_1.0-7              spatstat.utils_3.0-1
#  [35] DelayedArray_0.22.0       assertthat_0.2.1
#  [37] promises_1.2.0.1          scales_1.2.1
#  [39] ggraph_2.0.5              googlesheets4_1.0.0
#  [41] rgeos_0.5-9               beeswarm_0.4.0
#  [43] gtable_0.3.0              beachmat_2.12.0
#  [45] globals_0.15.1            goftest_1.2-3
#  [47] tidygraph_1.2.1           rlang_1.0.6
#  [49] splines_4.2.0             lazyeval_0.2.2
#  [51] gargle_1.2.0              spatstat.geom_3.0-3
#  [53] broom_1.0.0               reshape2_1.4.4
#  [55] abind_1.4-5               modelr_0.1.8
#  [57] backports_1.4.1           httpuv_1.6.5
#  [59] tools_4.2.0               ellipsis_0.3.2
#  [61] spatstat.core_2.4-4       RColorBrewer_1.1-3
#  [63] ggridges_0.5.3            Rcpp_1.0.10
#  [65] plyr_1.8.7                sparseMatrixStats_1.8.0
#  [67] zlibbioc_1.42.0           RCurl_1.98-1.10
#  [69] rpart_4.1.16              deldir_1.0-6
#  [71] pbapply_1.5-0             viridis_0.6.2
#  [73] cowplot_1.1.1             zoo_1.8-10
#  [75] haven_2.5.0               ggrepel_0.9.1
#  [77] cluster_2.1.3             fs_1.5.2
#  [79] scattermore_0.8           lmtest_0.9-40
#  [81] reprex_2.0.1              RANN_2.6.1
#  [83] googledrive_2.0.0         fitdistrplus_1.1-8
#  [85] hms_1.1.1                 mime_0.12
#  [87] xtable_1.8-4              readxl_1.4.0
#  [89] gridExtra_2.3             compiler_4.2.0
#  [91] KernSmooth_2.23-20        crayon_1.5.1
#  [93] htmltools_0.5.3           mgcv_1.8-40
#  [95] later_1.3.0               tzdb_0.3.0
#  [97] lubridate_1.8.0           DBI_1.1.3
#  [99] tweenr_1.0.2              dbplyr_2.2.1
# [101] MASS_7.3-56               Matrix_1.4-1
# [103] cli_3.4.1                 metapod_1.4.0
# [105] parallel_4.2.0            igraph_1.3.4
# [107] pkgconfig_2.0.3           plotly_4.10.0
# [109] spatstat.sparse_3.0-0     xml2_1.3.3
# [111] vipor_0.4.5               dqrng_0.3.0
# [113] XVector_0.36.0            rvest_1.0.2
# [115] digest_0.6.29             sctransform_0.3.3
# [117] RcppAnnoy_0.0.19          spatstat.data_3.0-0
# [119] cellranger_1.1.0          leiden_0.4.2
# [121] uwot_0.1.14               DelayedMatrixStats_1.18.0
# [123] shiny_1.7.2               gtools_3.9.4
# [125] nlme_3.1-157              lifecycle_1.0.3
# [127] jsonlite_1.8.0            BiocNeighbors_1.14.0
# [129] viridisLite_0.4.0         fansi_1.0.3
# [131] pillar_1.8.1              lattice_0.20-45
# [133] fastmap_1.1.0             httr_1.4.3
# [135] survival_3.3-1            glue_1.6.2
# [137] png_0.1-7                 bluster_1.6.0
# [139] ggforce_0.3.3             stringi_1.7.8
# [141] BiocSingular_1.12.0       irlba_2.3.5
# [143] future.apply_1.9.0
