#!/usr/bin/env Rscript
# File: cite_seq.r
# Author: Zhenhua Zhang
# E-mail: zhenhua.zhang217@gmail.com
# Created: Jun 06, 2023
# Updated: Mar 01, 2024
options(future.globals.maxSize=4096*1024^2)

suppressPackageStartupMessages({
  library(lobstr)
  library(tidyverse)
  library(data.table)
  library(ggsci)

  library(Seurat)
  library(Azimuth)
  library(SeuratData)
  library(CellMixS)
  library(SingleCellExperiment)
  library(scDblFinder)
})


proj_dir <- "~/Documents/projects/wp_pml"
deg_dir <- file.path(proj_dir, "outputs/analysis/CITE_seq/deg")
plot_dir <- file.path(proj_dir, "outputs/analysis/CITE_seq/plots")
object_dir <- file.path(proj_dir, "outputs/analysis/CITE_seq/objects")


# Cell type markers level 1 ()
# Shared markers
# TYROBP, (monocyte, NK cells), TRAC, (CD4+ T, CD8+ T), TMSB10, (CD4+ T, CD8+ T), CD3D, (CD4+ T, CD8+ T, Other T)
# CD3G, (CD4+ T, CD8+ T), NKG7, (NK, Other T), CST7, (NK, Other T), CD74, (DC, B), HLA-DQA1, (B, DC)
tar_cell_types <- c("Mono", "CD4 T", "CD8 T", "B", "NK", "DC")
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
## CITE-seq analysis
#
meta_data_path <- file.path(proj_dir, "misc/metadata/clinical_parameters.csv")
major_meta_tab <- fread(meta_data_path) %>%
  dplyr::select(Patient_ID, Age, Sex, BMI, Response) %>%
  dplyr::mutate(Patient_ID = stringr::str_remove(Patient_ID, "Cure-"))

# Sequencing meta data, samples per pool
other_meta_tab <- tibble::tribble(
  ~Patient_ID, ~SequencingPool, ~Timepoint, # ~Timepoint_raw,
  "PML0002", "CITEpool2", "3M", # "3 months",
  "PML0002", "CITEpool4", "BL", # "BL",
  "PML0009", "CITEpool3", "3M", # "3 months",
  "PML0009", "CITEpool2", "BL", # "BL",
  "PML0016", "CITEpool1", "BL", # "BL",
  "PML0035", "CITEpool3", "BL", # "BL",
  "PML0017", "CITEpool1", "6W", # "6 weeks",
  "PML0017", "CITEpool4", "BL", # "BL",
  "PML0022", "CITEpool1", "3M", # "3 months",
  "PML0022", "CITEpool2", "6W", # "6 weeks",
  "PML0022", "CITEpool3", "BL", # "BL",
  "PML0025", "CITEpool1", "6M", # "6 months",
  "PML0025", "CITEpool4", "BL", # "BL",
  "PML0008", "CITEpool4", "BL", # "BL",
  "PML0008", "CITEpool2", "6W", # "6 weeks",
  "PML0008", "CITEpool3", "3M", # "3 months",
)

# Load data and basic QC
selected_pools <- 1:4
pools <- paste0("CITEpool", selected_pools) %>% purrr::set_names(.)
bcmat_path <- file.path(proj_dir, "outputs/readcounts/CITE_seq", pools, "outs") %>% purrr::set_names(pools)
dmres_path <- file.path(proj_dir, "outputs/demultiplex/CITE_seq", pools, "vireo_outdir_ref/donor_ids.tsv") %>% purrr::set_names(pools)

pbmc_list <- lapply(pools, function(pn, .dmres_path, .bcmat_path, .overwrite) {
  pp <- file.path(.bcmat_path[[pn]], "filtered_feature_bc_matrix")
  obj_save_to <- file.path(object_dir, "raw", paste("pbmc", pn, "raw.rds", sep = "."))
  if (!file.exists(obj_save_to) || .overwrite) {
    cat("[I]: Creating from raw data ...", pn, "\n")
    dm_tab <- fread(.dmres_path[pn]) %>%
      dplyr::select(Cellbarcode = cell, Vireo_assignment = donor_id, Vireo_prob_max = prob_max, Vireo_prob_doublet = prob_doublet, Vireo_dbl_llr = doublet_logLikRatio) %>%
      dplyr::mutate(Vireo_assignment = dplyr::case_when(
        Vireo_assignment %in% c("doublet", "unassigned") ~ Vireo_assignment,
        TRUE ~ paste0("PML", stringr::str_remove(Vireo_assignment, "PML") %>% stringr::str_pad(4, pad = "0"))
      ))

    # Load CellRanger counts results and create RNA and ADT assaies
    pm <- Read10X(pp)
    po <- CreateSeuratObject(count = pm$`Gene Expression`, project = "PMLCITEseq")
    po[["ADT"]] <- CreateAssayObject(count = pm$`Antibody Capture`)

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
    p <- VlnPlot(object = po, features = c("nCount_RNA", "nFeature_RNA", "nCount_ADT", "nFeature_ADT", "percent_mt", "percent_rb"), ncol = 5, pt.size = 0)
    save_to <- file.path(plot_dir, "quality_control", paste("control_metrics", pn, "violin.pdf", sep = "."))
    ggsave(save_to, p, width = 12, height = 8)

    p <- FeatureScatter(po, feature1 = "nCount_RNA", feature2 = "percent_mt") +
      FeatureScatter(po, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
      FeatureScatter(po, feature1 = "nCount_ADT", feature2 = "nFeature_ADT")
    save_to <- file.path(plot_dir, "quality_control", paste("control_metrics", pn, "feature_scatter.pdf", sep = "."))
    ggsave(save_to, p, width = 12, height = 7)

    # Dump into disk
    cat("[I]: Pool", pn, "has", ncol(po), "cells and", nrow(po), "genes.\n")
    cat("[I]: Dumping into disk ...", pn, "\n"); saveRDS(po, obj_save_to)
  } else {
    cat("[I]: Loading from disk ...", pn, "\n"); po <- readRDS(obj_save_to)
  }

  invisible(po)
}, .dmres_path = dmres_path, .bcmat_path = bcmat_path, .overwrite = FALSE)


# Filtering cells using soft threshold quantile(probs = c(0.05, 0.95)) on each parameter.
soft_threshold <- tibble::tribble(
  ~SequencingPool, ~min_ncount_rna, ~max_ncount_rna, ~min_nfeature_rna, ~max_nfeature_rna, ~max_percent_mt,
  "CITEpool1", 2087, 17112, 986, 4190, 9.4,
  "CITEpool2", 1811, 11152, 773, 3380, 9.8,
  "CITEpool3", 1763, 20000, 574, 3742, 9.5,
  "CITEpool4", 2037, 18513, 1041, 4393, 9.4,
)

# selected_donors <- other_meta_tab$Patient_ID %>% unique(); save_token <- "all_samples"
selected_donors <- c("PML0002", "PML0009", "PML0017", "PML0022", "PML0025")

pbmc <- lapply(pools, function(pn, .pbmc, .threshold, .cc_genes) {
  # pn <- "CITEpool1"
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
    dplyr::filter(min_ncount_rna <= nCount_RNA, nCount_RNA <= max_ncount_rna, min_nfeature_rna <= nFeature_RNA, nFeature_RNA <= max_nfeature_rna, 100 <= nFeature_ADT, 500 <= nCount_ADT, nCount_ADT <= 20000, percent_mt < max_percent_mt, percent_rb < 50) %>%
    dplyr::filter(Vireo_assignment %in% selected_donors) %>%
    # dplyr::filter(!(Vireo_assignment == "PML0002" & Timepoint == "3M")) %>%
    rownames()

  tar_features <- c(rownames(po@assays$RNA), rownames(po@assays$ADT)) %>% purrr::discard(~stringr::str_detect(.x, "^MT-|^RP[LS]"))

  po <- po[tar_features, tar_cells] %>%
    NormalizeData(assay = "RNA", verbose = FALSE) %>%
    CellCycleScoring(assay = "RNA", s.features = .cc_genes$s.genes, g2m.features = .cc_genes$g2m.genes, set.ident = TRUE, verbose = FALSE) %>%
    SCTransform(assay = "RNA", vst.flavor = "v2", vars.to.regress = c("percent_mt", "S.Score", "G2M.Score"), method = "glmGamPoi", variable.features.n = 3000, verbose = FALSE) %>%
    NormalizeData(assay = "ADT", normalization.method = "CLR", margin = 2, verbose = FALSE) %>%
    RunPCA(assay = "SCT", verbose = FALSE)

  invisible(po)
}, .pbmc = pbmc_list, .threshold = soft_threshold, .cc_genes = cc.genes)

# Integrate multiple pools
features <- SelectIntegrationFeatures(pbmc, verbose = FALSE)
pbmc <- PrepSCTIntegration(pbmc, anchor.features = features, verbose = FALSE)
anchors <- FindIntegrationAnchors(pbmc, normalization.method = "SCT", anchor.features = features, reduction = "rpca", dims = 1:30, k.anchor = 20, verbose = FALSE)
pbmc_int <- IntegrateData(anchors, normalization.method = "SCT", dims = 1:30, verbose = FALSE)

# Do PCA and regress out cell cycle scores
pbmc_int <- RunPCA(pbmc_int, assay = "integrated", reduction.name = "intpca", verbose = FALSE) %>%
  ScaleData(assay = "ADT", vars.to.regress = c("percent_mt", "S.Score", "G2M.Score"), verbose = FALSE) %>%
  RunPCA(assay = "ADT", features = rownames(pbmc_int[["ADT"]]), reduction.name = "adtpca", verbose = FALSE)

# Find neighbors using multiple modal data and do UMAP
pbmc_int <- FindMultiModalNeighbors(pbmc_int, reduction.list = list("intpca", "adtpca"), dims.list = list(1:30, 1:30), modality.weight.name = "RNA.weight", verbose = FALSE)
pbmc_int <- RunUMAP(pbmc_int, nn.name = "weighted.nn", reduction = "intpca", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_", dims = 1:30, verbose = FALSE)
pbmc_int <- FindClusters(pbmc_int, graph.name = "wsnn", resolution = 0.1, verbose = FALSE)

# Identify doublets
DefaultAssay(pbmc_int) <- "SCT"
pbmc_cmb_sce <- as.SingleCellExperiment(pbmc_int)
reducedDim(pbmc_cmb_sce, "UMAP") <- pbmc_int@reductions$wnn.umap@cell.embeddings
reducedDim(pbmc_cmb_sce, "PCA") <- pbmc_int@reductions$intpca@cell.embeddings
pbmc_cmb_sce <- scDblFinder(pbmc_cmb_sce, clusters = TRUE)
pbmc_int$scDblFinder_class <- pbmc_cmb_sce$scDblFinder.class
pbmc_int$scDblFinder_score <- pbmc_cmb_sce$scDblFinder.score


#
## Result overview
#
# Estimate the entropy by CellMixS
entropy_plot_save_to <- file.path(plot_dir, "integrated", paste0("pbmc.cite_seq.integrated.umap_entropy.pdf"))
if (!file.exists(entropy_plot_save_to)) {
  pbmc_cmb_sce <- as.SingleCellExperiment(pbmc_int)
  reducedDim(pbmc_cmb_sce, "UMAP") <- pbmc_int@reductions$wnn.umap@cell.embeddings
  reducedDim(pbmc_cmb_sce, "PCA") <- pbmc_int@reductions$intpca@cell.embeddings
  pbmc_cmb_sce <- entropy(pbmc_cmb_sce, "SequencingPool", k = 20)
  pdf(entropy_plot_save_to, width = 12, height = 5)
  visOverview(pbmc_cmb_sce, "SequencingPool", dim_red = "UMAP", metric = "entropy")
  dev.off()
}

# Check the clusters, PCA
for (ppca in c("intpca", "adtpca")) {
  pca_plot_save_to <- file.path(plot_dir, "integrated", paste0("pbmc.cite_seq.integrated.", ppca, ".pdf"))
  p <- DimPlot(pbmc_int, reduction = ppca)
  ggsave(pca_plot_save_to, plot = p)
}

# Check the clusters, UMAP
for (per_group in c("seurat_clusters", "SequencingPool", "Response", "Sex", "Vireo_assignment", "Timepoint")) {
  nm <- stringr::str_replace_all(per_group, "\\.", "_")
  umap_save_to <- file.path(plot_dir, "integrated", paste0("pbmc.cite_seq.integrated.umap_by_", nm, ".pdf"))
  p <- DimPlot(pbmc_int, reduction = "wnn.umap", group.by = per_group)
  plot_wd <- ifelse(per_group == "predicted.celltype.l2", 10, 7)
  ggsave(umap_save_to, plot = p, width = plot_wd)
}

# Cell type markers
for (per_group in c("seurat_clusters")) {
  nm <- stringr::str_replace_all(per_group, "\\.", "_")
  save_to <- file.path(plot_dir, "integrated", paste0("pbmc.cite_seq.celltype_markers.dotplot_by_", nm, ".pdf"))
  p <- DotPlot(pbmc_int, features = celltype_markers, group.by = per_group, cluster.idents = TRUE) + coord_flip() + RotatedAxis()
  plot_wd <- pbmc_int[[per_group]] %>% unique() %>% nrow() %>% `*`(0.4) %>% `+`(2)
  ggsave(save_to, p, width = plot_wd, height = 12)
}

# Expression distribution of cell type marker genes
p <- FeaturePlot(pbmc_int, features = celltype_markers, order = TRUE, reduction = "wnn.umap", ncol = 5) &
  labs(x = NULL, y = NULL) &
  scale_y_continuous(expand = c(0.01, 0.01)) &
  scale_x_continuous(expand = c(0.01, 0.01)) &
  theme(legend.position = "none", axis.text = element_blank(), axis.ticks = element_blank(), axis.line = element_line(linetype = "dotted"), plot.title = element_text(size = 8))
save_to <- file.path(plot_dir, "integrated/pbmc.cite_seq.integrated.cell_type_markers.wnn_umap.pdf")
ggsave(save_to, plot = p, width = 10, height = 20)

# Check example
save_to <- file.path(plot_dir, "integrated", paste0("pbmc.cite_seq.integrated.CD16_gex.feature_plot.pdf"))
p1 <- FeaturePlot(pbmc_int, "Hu.CD16", cols = c("lightgrey", "darkgreen"), reduction = "wnn.umap") + ggtitle("CD16 protein")
p2 <- FeaturePlot(pbmc_int, "FCGR3A", reduction = "wnn.umap") + ggtitle("FCGR3A(CD16) expression")
p <- p1 | p2
ggsave(save_to, plot = p, width = 12, height = 6)

# Check and remove doublets
Idents(pbmc_int) <- "scDblFinder_class"
p <- DimPlot(pbmc_int, split.by = "scDblFinder_class", reduction = "wnn.umap", order = TRUE) & NoLegend()
save_to <- file.path(plot_dir, "integrated/pbmc.cite_seq.integrated.scdblfinder_class.umap.pdf")
ggsave(save_to, plot = p, width = 6, height = 4)

# Remove doublets
tar_cells <- colnames(pbmc_int)[pbmc_int$scDblFinder_class == "singlet"]
pbmc_int <- pbmc_int[, tar_cells]

# Prepare for DEG analysis
pbmc_int <- PrepSCTFindMarkers(pbmc_int, verbose = FALSE)

# Save the integrated object
obj_save_to <- file.path(object_dir, "integrated", paste0("pbmc.cite_seq.integrated.pca_umap_clustered.annotated.rds"))
saveRDS(pbmc_int, obj_save_to)
