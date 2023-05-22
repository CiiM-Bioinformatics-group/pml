#!/usr/bin/env Rscript

# TODO:
#   [ ] 1. Now many cells were annotated as B cells, test different filtering paramters


options(stringsAsFactors = FALSE, data.table.verbose = FALSE, error = rlang::last_trace)
suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(Seurat)
  library(Azimuth)
  library(SeuratData)
})


projdir <- "~/Documents/projects/wp_pml"
overwrite <- FALSE

# Load data
pools <- c("CITEpool1", "CITEpool2", "CITEpool3", "CITEpool4")
cite_seq_path <- file.path(projdir, "outputs/CITE_seq", pools, "outs/filtered_feature_bc_matrix")
names(cite_seq_path) <- pools

pbmc <- lapply(cite_seq_path, function(pp, .overwrite) {
  pn <- stringr::str_split(pp, pattern = "/", simplify = TRUE) %>% purrr::keep(~stringr::str_detect(.x, "CITEpool"))
  pm <- Read10X(pp)
  po <- CreateSeuratObject(count = pm$`Gene Expression`)
  po[["ADT"]] <- CreateAssayObject(count = pm$`Antibody Capture`, sep = c(":", "-"))

  save_to <- file.path(projdir, "outputs/analysis/cite_seq/objects", paste("pbmc", pn, "raw.rds", sep = "."))
  if (!file.exists(save_to) || .overwrite) saveRDS(po, save_to)

  DefaultAssay(po) <- "RNA"
  po[["percent.mt"]] <- PercentageFeatureSet(po, pattern = "^MT-")

  p <- VlnPlot(object = po, features = c("nCount_RNA", "nFeature_RNA", "nCount_ADT", "nFeature_ADT", "percent.mt"), ncol = 5, pt.size = 0)
  save_to <- file.path(projdir, "outputs/analysis/cite_seq/quality_control", paste("control_metrics", pn, "violin.pdf", sep = "."))
  ggsave(save_to, p, width = 12, height = 8)

  p <- FeatureScatter(po, feature1 = "nCount_RNA", feature2 = "percent.mt") +
    FeatureScatter(po, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
    FeatureScatter(po, feature1 = "nCount_ADT", feature2 = "nFeature_ADT")
  save_to <- file.path(projdir, "outputs/analysis/cite_seq/quality_control", paste("control_metrics", pn, "feature_scatter.pdf", sep = "."))
  ggsave(save_to, p, width = 12, height = 7)

  po
}, .overwrite = FALSE)



hard_threshold <- tibble::tribble(
  ~pool, ~min_ncount_rna, ~max_ncount_rna, ~min_nfeature_rna, ~max_nfeature_rna, ~max_percent_mt,
  "CITEpool1", 500, 20000, 1000, 5000, 7.5,
  "CITEpool2", 500, 25000, 750, 5000, 7.5,
  "CITEpool3", 500, 20000, 750, 4500, 7.5,
  "CITEpool4", 500, 25000, 750, 5250, 7.5,
)


pools <- names(pbmc)
names(pools) <- pools
pbmc <- pools %>%
  lapply(function(po, .pbmc, .hard_threshold) {
    params <- .hard_threshold %>% dplyr::filter(pool == po)
    min_ncount_rna <- params["min_ncount_rna"] %>% as.integer()
    max_ncount_rna <- params["max_ncount_rna"] %>% as.integer()
    min_nfeature_rna <- params["min_nfeature_rna"] %>% as.integer()
    max_nfeature_rna <- params["max_nfeature_rna"] %>% as.integer()
    max_percent_mt <- params["max_percent_mt"] %>% as.double()

    tar_cells <- .pbmc[[po]]@meta.data %>%
      as.data.frame() %>%
      dplyr::filter(
        min_ncount_rna <= nCount_RNA, nCount_RNA <= max_ncount_rna,
        min_nfeature_rna <= nFeature_RNA, nFeature_RNA <= max_nfeature_rna,
        100 <= nFeature_ADT, 500 <= nCount_ADT, nCount_ADT <= 20000,
        percent.mt < max_percent_mt
      ) %>%
      rownames()

    .pbmc[[po]]@meta.data$pool <- po
    NormalizeData(.pbmc[[po]][, tar_cells], assay = "RNA") %>%
      NormalizeData(normalization.method = "CLR", margin = 2, assay = "ADT")
}, .pbmc = pbmc, .hard_threshold = hard_threshold)


#
## Integrate multiple pools by features from RNA assay.
#
pbmc <- lapply(pbmc, function(po) { DefaultAssay(po) <- "RNA"; po})
features <- SelectIntegrationFeatures(pbmc)
anchors <- FindIntegrationAnchors(pbmc, anchor.features = features)
pbmc_cmb <- IntegrateData(anchorset = anchors)

save_to <- file.path(projdir, "outputs/analysis/cite_seq/objects/pbmc.qc.integrated.rds")
if (!file.exists(save_to) || overwrite) {
  saveRDS(pbmc_cmb, save_to)
} else {
  pbmc_cmb <- readRDS(save_to)
}

# Do PCA
DefaultAssay(pbmc_cmb) <- "integrated"
pbmc_cmb <- ScaleData(pbmc_cmb)
pbmc_cmb <- RunPCA(pbmc_cmb, npcs = 30, reduction.name = "int_pca")

DefaultAssay(pbmc_cmb) <- "ADT"
VariableFeatures(pbmc_cmb) <- rownames(pbmc_cmb[["ADT"]])
pbmc_cmb <- ScaleData(pbmc_cmb)
pbmc_cmb <- RunPCA(pbmc_cmb, npcs = 30, reduction.name = "adt_pca")

save_to <- file.path(projdir, "outputs/analysis/cite_seq/objects/pbmc.qc.integrated.pca.rds")
if (!file.exists(save_to) || overwrite) {
  saveRDS(pbmc_cmb, save_to)
} else {
  pbmc_cmb <- readRDS(save_to)
}

p <- DimPlot(pbmc_cmb, reduction = "int_pca")
save_to <- file.path(projdir, "outputs/analysis/cite_seq/overview/pbmc.qc.integrated.int_pca.pdf")
ggsave(save_to, plot = p)

p <- DimPlot(pbmc_cmb, reduction = "adt_pca")
save_to <- file.path(projdir, "outputs/analysis/cite_seq/overview/pbmc.qc.integrated.adt_pca.pdf")
ggsave(save_to, plot = p)

# Find neighbors and do UMAP
pbmc_cmb <- FindMultiModalNeighbors(pbmc_cmb, reduction.list = list("int_pca", "adt_pca"), dims.list = list(1:30, 1:30), modality.weight.name = "RNA.weight")
pbmc_cmb <- RunUMAP(pbmc_cmb, nn.name = "weighted.nn", reduction = "int_pca", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_", dims = 1:30)

save_to <- file.path(projdir, "outputs/analysis/cite_seq/objects/pbmc.qc.integrated.pca.clustered.rds")
if (!file.exists(save_to) || overwrite) {
  saveRDS(pbmc_cmb, save_to)
} else {
  pbmc_cmb <- readRDS(save_to)
}

# Create clusters
# pbmc_cmb <- FindNeighbors(pbmc_cmb, reduction = "pca", dims = 1:30)
pbmc_cmb <- FindClusters(pbmc_cmb, graph.name = "wsnn", resolution = 0.1)

Idents(pbmc_cmb) <- pbmc_cmb$seurat_clusters
p <- DimPlot(pbmc_cmb, reduction = "wnn.umap")
save_to <- file.path(projdir, "outputs/analysis/cite_seq/overview/pbmc.qc.integrated.pca.clustered.umap_by_cluster.pdf")
ggsave(save_to, plot = p)

Idents(pbmc_cmb) <- pbmc_cmb$pool
p <- DimPlot(pbmc_cmb, reduction = "wnn.umap")
save_to <- file.path(projdir, "outputs/analysis/cite_seq/overview/pbmc.qc.integrated.pca.clustered.umap_by_pool.pdf")
ggsave(save_to, plot = p)


# Check example
DefaultAssay(pbmc_cmb) <- "ADT"
p1 <- FeaturePlot(pbmc_cmb, "Hu.CD16", cols = c("lightgrey", "darkgreen"), reduction = "wnn.umap") + ggtitle("CD16 protein")

DefaultAssay(pbmc_cmb) <- "RNA"
p2 <- FeaturePlot(pbmc_cmb, "FCGR3A", reduction = "wnn.umap") + ggtitle("FCGR3A(CD16) expression")

p <- p1 | p2
save_to <- file.path(projdir, "outputs/analysis/cite_seq/overview/pbmc.integrated.CD16.gex.feature_plot.pdf")
ggsave(save_to, plot = p, width = 10, height = 6)


# Annotate the cells by Azimuth
pbmc_cmb <- RunAzimuth(pbmc_cmb, reference = "pbmcref")

save_to <- file.path(projdir, "outputs/analysis/cite_seq/objects/pbmc.qc.integrated.pca.clustered.annotated.rds")
if (!file.exists(save_to) || overwrite) {
  saveRDS(pbmc_cmb, save_to)
} else {
  pbmc_cmb <- readRDS(save_to)
}

Idents(pbmc_cmb) <- pbmc_cmb$predicted.celltype.l1
p <- DimPlot(pbmc_cmb, reduction = "wnn.umap", label = TRUE) + NoLegend()
save_to <- file.path(projdir, "outputs/analysis/cite_seq/overview/pbmc.qc.integrated.pca.clustered.umap_by_celltype_l1.pdf")
ggsave(save_to, plot = p)

p <- FeaturePlot(pbmc_cmb, features = c("CD79A", "CD79B", "MS4A1", "CD74", "BANK1", "RALGPS2"), ncol = 3) & NoAxes()
save_to <- file.path(projdir, "outputs/analysis/cite_seq/overview/pbmc.qc.integrated.pca.clustered.b_cell.feature_plot.pdf")
ggsave(save_to, plot = p, width = 12)


p <- FeaturePlot(pbmc_cmb, features = c("IL7R", "MAL", "LTB", "CD4", "LDHB", "TPT1"), ncol = 3) & NoAxes()
save_to <- file.path(projdir, "outputs/analysis/cite_seq/overview/pbmc.qc.integrated.pca.clustered.cd4_t.feature_plot.pdf")
ggsave(save_to, plot = p, width = 12)

p <- FeaturePlot(pbmc_cmb, features = c("CD8B", "CD8A", "CD3D", "TMSB10", "HCST", "CD3G"), ncol = 3) & NoAxes()
save_to <- file.path(projdir, "outputs/analysis/cite_seq/overview/pbmc.qc.integrated.pca.clustered.cd8_t.feature_plot.pdf")
ggsave(save_to, plot = p, width = 12)

p <- FeaturePlot(pbmc_cmb, features = c("CTSS", "FCN1", "NEAT1", "LYZ", "PSAP", "S100A9"), ncol = 3) & NoAxes()
save_to <- file.path(projdir, "outputs/analysis/cite_seq/overview/pbmc.qc.integrated.pca.clustered.monocyte.feature_plot.pdf")
ggsave(save_to, plot = p, width = 12)

# Add demultiplex information
new_meta_tab <- fread()



#
## Merge multiple pools, does not integrate different pool properly.
#
if (FALSE) {
  pbmc_merge <- merge(
    pbmc$CITEpool1, y = c(pbmc$CITEpool2, pbmc$CITEpool3, pbmc$CITEpool4),
    add.cell.ids = c("pool1", "pool2", "pool3", "pool4"), project = "PML_CITEseq"
  )

  DefaultAssay(pbmc_merge) <- "RNA"
  pbmc_merge <- ScaleData(pbmc_merge)
  pbmc_merge <- FindVariableFeatures(pbmc_merge)
  pbmc_merge <- RunPCA(pbmc_merge, npcs = 30, reduction.name = "int_pca")

  DefaultAssay(pbmc_merge) <- "ADT"
  VariableFeatures(pbmc_merge) <- rownames(pbmc_merge[["ADT"]])
  pbmc_merge <- ScaleData(pbmc_merge)
  pbmc_merge <- RunPCA(pbmc_merge, npcs = 30, reduction.name = "adt_pca")

  save_to <- file.path(projdir, "outputs/analysis/cite_seq/objects/pbmc.qc.merged.pca.rds")
  if (!file.exists(save_to) || overwrite) {
    saveRDS(pbmc_merge, save_to)
  } else {
    pbmc_merge <- readRDS(save_to)
  }

  p <- DimPlot(pbmc_merge, reduction = "int_pca")
  save_to <- file.path(projdir, "outputs/analysis/cite_seq/overview/pbmc.qc.merged.int_pca.pdf")
  ggsave(save_to, plot = p)

  p <- DimPlot(pbmc_merge, reduction = "adt_pca")
  save_to <- file.path(projdir, "outputs/analysis/cite_seq/overview/pbmc.qc.merged.adt_pca.pdf")
  ggsave(save_to, plot = p)

  # Find neighbors and do UMAP
  pbmc_merge <- FindMultiModalNeighbors(pbmc_merge, reduction.list = list("int_pca", "adt_pca"), dims.list = list(1:30, 1:30), modality.weight.name = "RNA.weight")
  pbmc_merge <- RunUMAP(pbmc_merge, nn.name = "weighted.nn", reduction = "int_pca", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_", dims = 1:30)

  save_to <- file.path(projdir, "outputs/analysis/cite_seq/objects/pbmc.qc.merged.pca.clustered.rds")
  if (!file.exists(save_to) || overwrite) {
    saveRDS(pbmc_merge, save_to)
  } else {
    pbmc_merge <- readRDS(save_to)
  }

  # Create clusters
  # pbmc_merge <- FindNeighbors(pbmc_merge, reduction = "pca", dims = 1:30)
  pbmc_merge <- FindClusters(pbmc_merge, graph.name = "wsnn", resolution = 0.1)

  Idents(pbmc_merge) <- pbmc_merge$seurat_clusters
  p <- DimPlot(pbmc_merge, reduction = "wnn.umap")
  save_to <- file.path(projdir, "outputs/analysis/cite_seq/overview/pbmc.qc.merged.pca.clustered.umap_by_cluster.pdf")
  ggsave(save_to, plot = p)

  Idents(pbmc_merge) <- pbmc_merge$pool
  p <- DimPlot(pbmc_merge, reduction = "wnn.umap")
  save_to <- file.path(projdir, "outputs/analysis/cite_seq/overview/pbmc.qc.merged.pca.clustered.umap_by_pool.pdf")
  ggsave(save_to, plot = p)


  # Check example
  DefaultAssay(pbmc_merge) <- "ADT"
  p1 <- FeaturePlot(pbmc_merge, "Hu.CD16", cols = c("lightgrey", "darkgreen"), reduction = "wnn.umap") + ggtitle("CD16 protein")

  DefaultAssay(pbmc_merge) <- "RNA"
  p2 <- FeaturePlot(pbmc_merge, "FCGR3A", reduction = "wnn.umap") + ggtitle("FCGR3A(CD16) expression")

  p <- p1 | p2
  save_to <- file.path(projdir, "outputs/analysis/cite_seq/overview/pbmc.merged.CD16.gex.feature_plot.pdf")
  ggsave(save_to, plot = p, width = 10, height = 6)


  # annotate the cells by Azimuth
  pbmc_merge <- RunAzimuth(pbmc_merge, reference = "pbmcref")

  save_to <- file.path(projdir, "outputs/analysis/cite_seq/objects/pbmc.qc.merged.pca.clustered.annotated.rds")
  if (!file.exists(save_to) || overwrite) {
    saveRDS(pbmc_merge, save_to)
  } else {
    pbmc_merge <- readRDS(save_to)
  }

  Idents(pbmc_merge) <- pbmc_merge$predicted.celltype.l1
  p <- DimPlot(pbmc_merge, reduction = "wnn.umap", label = TRUE) + NoLegend()
  save_to <- file.path(projdir, "outputs/analysis/cite_seq/overview/pbmc.qc.merged.pca.clustered.umap_by_celltype_l1.pdf")
  ggsave(save_to, plot = p)

}



#
## Aggregated results by cellranger, does not correct pool batch properly.
#
if (FALSE) {
  cite_seq_path <- file.path(projdir, "outputs/CITE_seq/All/outs/count/filtered_feature_bc_matrix")
  pbmc_attr_mat <- Read10X(cite_seq_path)
  pbmc_aggr <- CreateSeuratObject(count = pbmc_attr_mat$`Gene Expression`)
  pbmc_aggr[["ADT"]] <- CreateAssayObject(count = pbmc_attr_mat$`Antibody Capture`, sep = c(":", "-"))

  save_to <- file.path(projdir, "outputs/analysis/cite_seq/objects/pbmc.aggr.raw.rds")
  if (!file.exists(save_to) || overwrite) saveRDS(pbmc_aggr, save_to)

  DefaultAssay(pbmc_aggr) <- "RNA"
  pbmc_aggr[["percent.mt"]] <- PercentageFeatureSet(pbmc_aggr, pattern = "^MT-")

  pbmc_aggr <- subset(x = pbmc_aggr, subset = 5e2 < nCount_RNA & nCount_RNA < 2e4 & 5e2 <= nFeature_RNA & nFeature_RNA <= 6e3 & percent.mt < 7.5)
  save_to <- file.path(projdir, "outputs/analysis/cite_seq/pbmc.aggr.filtered.rds")
  if (!file.exists(save_to) || overwrite) {
    saveRDS(pbmc_aggr, save_to)
  } else {
    pbmc_aggr <- readRDS(save_to)
    DefaultAssay(pbmc_aggr)
  }

  DefaultAssay(pbmc_aggr) <- "RNA"
  pbmc_aggr <- NormalizeData(pbmc_aggr)
  pbmc_aggr <- FindVariableFeatures(pbmc_aggr)
  top10 <- head(VariableFeatures(pbmc_aggr), 10)

  p1 <- VariableFeaturePlot(pbmc_aggr)
  p2 <- LabelPoints(plot = p1, points = top10, repel = TRUE)
  p <- p1 + p2
  save_to <- file.path(projdir, "outputs/analysis/cite_seq/overview/aggr.variable_features_top10.pdf")
  ggsave(save_to, plot = p, width = 12)

  pbmc_aggr <- ScaleData(pbmc_aggr)

  pbmc_aggr <- RunPCA(pbmc_aggr, verbose = FALSE)
  p <- DimPlot(pbmc_aggr, reduction = "pca")
  save_to <- file.path(projdir, "outputs/analysis/cite_seq/overview/aggr.pca.gex.pdf")
  ggsave(save_to, plot = p)

  pbmc_aggr <- FindNeighbors(pbmc_aggr, dims = 1:30)
  pbmc_aggr <- FindClusters(pbmc_aggr, resolution = .5, verbose = FALSE)
  pbmc_aggr@meta.data$pool <- pbmc_aggr@meta.data %>% rownames() %>% stringr::str_extract("[0-9]$")

  pbmc_aggr <- RunUMAP(pbmc_aggr, dims = 1:30)
  p <- DimPlot(pbmc_aggr, reduction = "umap", group.by = "pool")
  save_to <- file.path(projdir, "outputs/analysis/cite_seq/overview/aggr.umap.gex.pdf")
  ggsave(save_to, plot = p)

  save_to <- file.path(projdir, "outputs/analysis/cite_seq/pbmc.aggr.filtered.pca.umap.clustered.rds")
  if (!file.exists(save_to) || overwrite) {
    saveRDS(pbmc_aggr, save_to)
  } else {
    pbmc_aggr <- readRDS(save_to)
    DefaultAssay(pbmc_aggr)
  }
}
