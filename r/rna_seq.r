#!/usr/bin/env Rscript
# File: rna_seq.r
# Author: Zhenhua Zhang
# E-mail: zhenhua.zhang217@gmail.com
# Created: Aug 21, 2023
# Updated: Mar 01, 2024

suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(ggsci)

  library(Seurat)
  library(Azimuth)
  library(SeuratData)
  library(CellMixS)
  library(SingleCellExperiment)
  library(scDblFinder)
})


projdir <- "~/Documents/projects/wp_pml"
deg_dir <- file.path(projdir, "outputs/analysis/RNA_seq/deg")
plot_dir <- file.path(projdir, "outputs/analysis/RNA_seq/plots")
object_dir <- file.path(projdir, "outputs/analysis/RNA_seq/objects")


# Cell type markers level 1 ()
# Shared markers
# TYROBP, (monocyte, NK cells), TRAC, (CD4+ T, CD8+ T), TMSB10, (CD4+ T, CD8+ T), CD3D, (CD4+ T, CD8+ T, Other T)
# CD3G, (CD4+ T, CD8+ T), NKG7, (NK, Other T), CST7, (NK, Other T), CD74, (DC, B), HLA-DQA1, (B, DC)
celltype_markers <- c(
  "CTSS", "FCN1", "LYZ", "PSAP", "S100A9", "AIF1", "SERPINA1", "CD14", "FCGR3A", # "TYROBP", "NEAT1", "MNDA", Monocytes
  "IL7R", "MAL", "LTB", "CD4", "LDHB", "CD3D", "CD3G", "TRAC", # "TMSB10", "TPT1", CD4+ T
  "CD8B", "CD8A", "HCST", "LINC02446", "CTSW", "CD3E", # "TMSB10", "CD3G", "CD3D", "TRAC", CD8+ T
  "TRDC", "GZMK", "KLRB1", "TRGC2", "LYAR", "KLRG1", "GZMA", # "CD3D", "CST7", "NKG7", Other T
  "KLRD1", "GNLY", "PRF1", "CD247", "KLRF1", "GZMB", # "TYROBP", "FCER1G", "CST7", NK
  "CD79A", "CD79B", "RALGPS2", "MS4A1", "BANK1", "TNFRSF13C", "IGHM", "MEF2C", # "CD74", "HLA-DQA1", B
  "HLA-DPA1", "HLA-DPB1", "CCDC88A", "HLA-DRA", "HLA-DMA", "CST3", "HLA-DQB1", "HLA-DRB1" # "CD74", "HLA-DQA1", DC
)


#
## RNA-seq analysis
#
meta_data_path <- file.path(projdir, "misc/metadata/clinical_parameters.csv")
major_meta_tab <- fread(meta_data_path) %>%
  dplyr::select(Patient_ID, Age, Sex, BMI, Response) %>%
  dplyr::mutate(Patient_ID = stringr::str_remove(Patient_ID, "Cure-"))

# Sequencing meta data, samples per pool
other_meta_tab <- tibble::tribble(
  ~Patient_ID, ~SequencingPool, ~Timepoint, # ~Timepoint_raw,
  "PML0004", "RNAB2P1", "BL",
  "PML0063", "RNAB2P4", "BL",
  "PML0066", "RNAB2P2", "BL",
  "PML0070", "RNAB2P3", "BL",
  "PML0007", "RNAB2P3", "BL",
  "PML0013", "RNAB2P4", "BL",
  "PML0033", "RNAB2P3", "BL",
  "PML0045", "RNAB2P1", "BL",
  "PML0055", "RNAB2P2", "BL",
  "PML0058", "RNAB2P2", "BL",
  "PML0060", "RNAB2P4", "BL",
  "PML0061", "RNAB2P1", "BL",
)

# Load data and basic QC
selected_pools <- 1:4
pools <- paste0("RNAB2P", selected_pools) %>% rlang::set_names(.)
bcmat_path <- file.path(projdir, "outputs/readcounts/RNA_seq", pools, "outs") %>% rlang::set_names(pools)
dmres_path <- file.path(projdir, "outputs/demultiplex/RNA_seq", pools, "vireo_outdir_ref/donor_ids.tsv") %>% purrr::set_names(pools)

pbmc_list <- lapply(pools, function(pn, .dmres_path, .bcmat_path, .overwrite) {
  obj_save_to <- file.path(object_dir, "raw", paste0("pbmc.", pn, ".raw.rds"))
  pp <- file.path(.bcmat_path[pn], "filtered_feature_bc_matrix")
  if (!file.exists(obj_save_to) || .overwrite) {
    cat("[I]: Creating from raw data ...\n")
    dm_tab <- fread(.dmres_path[pn]) %>%
      dplyr::select(Cellbarcode = cell, Vireo_assignment = donor_id, Vireo_prob_max = prob_max, Vireo_prob_doublet = prob_doublet, Vireo_dbl_llr = doublet_logLikRatio)

    # Load CellRanger counts results and create RNA
    pm <- Read10X(pp)
    po <- CreateSeuratObject(pm, project = "PMLRNA", min.cells = 3, min.features = 200)

    # MT and RB gene ratios
    DefaultAssay(po) <- "RNA"
    po[["percent_mt"]] <- PercentageFeatureSet(po, pattern = "^MT-")
    po[["percent_rb"]] <- PercentageFeatureSet(po, pattern = "^RP[LS]")

    # Update meta data
    raw_barcodes <- colnames(po)
    extra_meta_tab <- read.csv(file.path(pp, "barcodes.tsv"), col.names = "Cellbarcode") %>%
      dplyr::mutate(SequencingPool = pn) %>%
      dplyr::left_join(dm_tab, by = "Cellbarcode") %>%
      dplyr::left_join(major_meta_tab, by = c("Vireo_assignment" = "Patient_ID")) %>%
      dplyr::left_join(other_meta_tab, by = c("Vireo_assignment" = "Patient_ID", "SequencingPool"))
    po@meta.data <- po@meta.data %>%
      dplyr::mutate(Cellbarcode = rownames(.) %>% stringr::str_remove("_[0-9]$"), SequencingPool = pn) %>%
      dplyr::left_join(extra_meta_tab, by = c("SequencingPool", "Cellbarcode")) %>%
      (function(tab) {rownames(tab) <- raw_barcodes; tab})

    # QC plots
    p <- VlnPlot(object = po, features = c("nCount_RNA", "nFeature_RNA", "percent_mt", "percent_rb"), ncol = 5, pt.size = 0)
    save_to <- file.path(plot_dir, "quality_control", paste0("control_metrics.", pn, ".violin.pdf"))
    ggsave(save_to, p, width = 12, height = 8)

    p <- FeatureScatter(po, feature1 = "nCount_RNA", feature2 = "percent_mt") +
      FeatureScatter(po, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

    save_to <- file.path(plot_dir, "quality_control", paste0("control_metrics.", pn, ".feature_scatter.pdf"))
    ggsave(save_to, p, width = 10, height = 7)

    cat("[I]: Pool", pn, "has", ncol(po), "cells and", nrow(po), "genes.\n")
    cat("[I]: Dumping into disk ...\n"); saveRDS(po, obj_save_to)
  } else {
    cat("[I]: Loading from disk ...\n"); po <- readRDS(obj_save_to)
  }

  invisible(po)
}, .dmres_path = dmres_path, .bcmat_path = bcmat_path, .overwrite = TRUE)


# Filtering cells using soft threshold quantile(probs = c(0.05, 0.95)) on each parameter.
soft_threshold <- tibble::tribble(
  ~SequencingPool, ~min_ncount_rna, ~max_ncount_rna, ~min_nfeature_rna, ~max_nfeature_rna, ~max_percent_mt,
  "RNAB2P1", 1378, 18271, 851, 4153, 9,
  "RNAB2P2", 550, 13984, 438, 3801, 9,
  "RNAB2P3", 1050, 28928, 383, 5248, 9,
  "RNAB2P4", 2084, 32602, 1026, 5534, 9,
)

# "PML0007", atypical B-cell proprotion, Responders
# "PML0045", atypical monocyte proprotion, Responders
# "PML0004", atypical B-cell proprotion, Non-responders
# "PML0066", atypical monocyte proprotion, Non-Responders
selected_donors <- c(
  "PML0063", "PML0070", "PML0004", "PML0066", # Non-responders
  "PML0007", "PML0045", "PML0013", "PML0033", "PML0055", "PML0058", "PML0060", "PML0061" # Responders
)
pbmc <- lapply(pools, function(pn, .pbmc, .threshold, .cc_genes) {
  params <- .threshold %>% dplyr::filter(SequencingPool == pn)
  min_ncount_rna <- params["min_ncount_rna"] %>% as.integer()
  max_ncount_rna <- params["max_ncount_rna"] %>% as.integer()
  min_nfeature_rna <- params["min_nfeature_rna"] %>% as.integer()
  max_nfeature_rna <- params["max_nfeature_rna"] %>% as.integer()
  max_percent_mt <- params["max_percent_mt"] %>% as.double()

  po <- .pbmc[[pn]]
  po$SequencingPool <- pn
  tar_cells <- po@meta.data %>%
    as.data.frame() %>%
    dplyr::filter(min_ncount_rna <= nCount_RNA, nCount_RNA <= max_ncount_rna, min_nfeature_rna <= nFeature_RNA, nFeature_RNA <= max_nfeature_rna, percent_mt < max_percent_mt, percent_rb < 50) %>%
    dplyr::filter(Vireo_assignment %in% selected_donors) %>%
    rownames()

  tar_features <- c(rownames(po@assays$RNA)) %>% purrr::discard(~stringr::str_detect(.x, "^MT-"))

  po <- po[tar_features, tar_cells] %>%
    NormalizeData(assay = "RNA", verbose = FALSE) %>%
    CellCycleScoring(assay = "RNA", s.features = .cc_genes$s.genes, g2m.features = .cc_genes$g2m.genes, set.ident = TRUE, verbose = FALSE) %>%
    SCTransform(assay = "RNA", vst.flavor = "v2", vars.to.regress = c("percent_mt", "S.Score", "G2M.Score"), method = "glmGamPoi", variable.features.n = 3000, verbose = FALSE) %>%
    RunPCA(assay = "SCT", verbose = FALSE)

  invisible(po)
}, .pbmc = pbmc_list, .threshold = soft_threshold, .cc_genes = cc.genes)

# Integrate multiple pools
features <- SelectIntegrationFeatures(pbmc)
pbmc <- PrepSCTIntegration(pbmc, anchor.features = features, verbose = FALSE)
anchors <- FindIntegrationAnchors(pbmc, normalization.method = "SCT", anchor.features = features, reduction = "rpca", dims = 1:30, k.anchor = 20, verbose = FALSE)
pbmc_int <- IntegrateData(anchors, normalization.method = "SCT", dims = 1:30, verbose = FALSE)

# Do PCA
pbmc_int <- RunPCA(pbmc_int, assay = "integrated", verbose = FALSE)

# Find neighbors using multiple modal data and do UMAP
pbmc_int <- RunUMAP(pbmc_int, assay = "integrated", dims = 1:30, verbose = FALSE)
pbmc_int <- FindNeighbors(pbmc_int, assay = "integrated", dims = 1:30, verbose = FALSE)
pbmc_int <- FindClusters(pbmc_int, graph.name = "integrated_nn", resolution = 0.1, verbose = FALSE)

# Identify doublets
DefaultAssay(pbmc_int) <- "SCT"
pbmc_cmb_sce <- as.SingleCellExperiment(pbmc_int)
reducedDim(pbmc_cmb_sce, "UMAP") <- pbmc_int@reductions$umap@cell.embeddings
reducedDim(pbmc_cmb_sce, "PCA") <- pbmc_int@reductions$pca@cell.embeddings
pbmc_cmb_sce <- scDblFinder(pbmc_cmb_sce, clusters = TRUE)
pbmc_int$scDblFinder_class <- pbmc_cmb_sce$scDblFinder.class
pbmc_int$scDblFinder_score <- pbmc_cmb_sce$scDblFinder.score

#
## Result overview
#
# Estimate the entropy by CellMixS
pbmc_cmb_sce <- entropy(pbmc_cmb_sce, "SequencingPool", k = 20)
save_to <- file.path(plot_dir, "integrated", "pbmc.rna_seq.integrated.umap_entropy.pdf")
pdf(save_to, width = 12, height = 5)
visOverview(pbmc_cmb_sce, group = "SequencingPool", metric = "entropy", dim_red = "UMAP")
dev.off()

# Check the clusters, PCA
p <- DimPlot(pbmc_int, reduction = "pca")
save_to <- file.path(plot_dir, "integrated/pbmc.rna_seq.integrated.pca.pdf")
ggsave(save_to, plot = p, width = 8, height = 8)

# Check the clusters, UMAP
for (per_group in c("seurat_clusters", "SequencingPool", "Response", "Sex", "Vireo_assignment", "Timepoint")) {
  p <- DimPlot(pbmc_int, reduction = "umap", group.by = per_group, pt.size = 0.75)
  nm <- stringr::str_replace_all(per_group, "\\.", "_")
  save_to <- file.path(plot_dir, "integrated", paste0("pbmc.rna_seq.integrated.umap_by_", nm, ".pdf"))
  plot_wd <- ifelse(per_group == "predicted.celltype.l2", 9, 6)
  ggsave(save_to, plot = p, width = plot_wd, height = 6)
}

# Cell type markers
for (per_group in c("seurat_clusters")) {
  p <- DotPlot(pbmc_int, features = celltype_markers, group.by = per_group, cluster.idents = TRUE) + coord_flip() + RotatedAxis()
  nm <- stringr::str_replace_all(per_group, "\\.", "_")
  save_to <- file.path(plot_dir, "integrated", paste0("pbmc.rna_seq.celltype_markers.dotplot_by_", nm, ".pdf"))
  plot_wd <- pbmc_int[[per_group]] %>% unique() %>% nrow() %>% `*`(0.4) %>% `+`(2)
  ggsave(save_to, p, width = plot_wd, height = 12)
}

# Expression distribution of cell type marker genes
p <- FeaturePlot(pbmc_int, features = celltype_markers, order = TRUE, reduction = "umap", ncol = 5) &
  labs(x = NULL, y = NULL) &
  scale_y_continuous(expand = c(0.01, 0.01)) &
  scale_x_continuous(expand = c(0.01, 0.01)) &
  theme(legend.position = "none", axis.text = element_blank(), axis.ticks = element_blank(), axis.line = element_line(linetype = "dotted"), plot.title = element_text(size = 8))
save_to <- file.path(plot_dir, "integrated/pbmc.rna_seq.integrated.cell_type_markers.umap.pdf")
ggsave(save_to, plot = p, width = 10, height = 20)

# Check and remove doublets
Idents(pbmc_int) <- "scDblFinder_class"
p <- DimPlot(pbmc_int, split.by = "scDblFinder_class", reduction = "umap", order = TRUE) & NoLegend()
save_to <- file.path(plot_dir, "integrated/pbmc.rna_seq.integrated.scdblfinder_class.umap.pdf")
ggsave(save_to, plot = p, width = 6, height = 4)

# Remove doublets
tar_cells <- colnames(pbmc_int)[pbmc_int$scDblFinder_class == "singlet"]
pbmc_int <- pbmc_int[, tar_cells]

# Prepare for DEG analysis
pbmc_int <- PrepSCTFindMarkers(pbmc_int, verbose = FALSE)

# Save the integrated object
save_to <- file.path(object_dir, "integrated/pbmc.rna_seq.integrated.pca_umap_clustered.annotated.rds")
saveRDS(pbmc_int, save_to)
