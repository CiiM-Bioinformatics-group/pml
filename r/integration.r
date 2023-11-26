#!/usr/bin/env Rscript
# File: integration_rna_atac.r
# Author: Zhenhua Zhang
# E-mail: zhenhua.zhang217@gmail.com
# Created: Aug 21, 2023
# Updated: Nov 22, 2023

suppressPackageStartupMessages({
  library(lobstr)

  library(data.table)
  library(tidyverse)
  library(ggsci)

  library(SingleCellExperiment)
  library(SeuratData)
  library(SeuratDisk)

  library(GenomicRanges)
  library(EnsDb.Hsapiens.v86)
  library(BSgenome.Hsapiens.UCSC.hg38)

  library(Seurat)
  library(Signac)
})


proj_dir <- "~/Documents/projects/wp_pml"
plot_dir <- file.path(proj_dir, "outputs/analysis/integrated/plots")
object_dir <- file.path(proj_dir, "outputs/analysis/integrated/objects")


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


# Load GEX reference dataset from Hao etal 2020
pbmc_reference <- file.path(proj_dir, "inputs/references/Hao_etal_2020/pbmc_multimodal.h5seurat") %>% LoadH5Seurat()
Misc(pbmc_reference[["wnn.umap"]], slot = "model")$num_precomputed_nns <- 1
for (per_idents in c(paste0("celltype.", c("l1", "l2")))) {
  Idents(pbmc_reference) <- per_idents
  p <- DimPlot(pbmc_reference, reduction = "wnn.umap", label = TRUE, repel = TRUE) + NoLegend()
  file.path(plot_dir, paste0("pbmc_reference_rnaseq.umap.", per_idents, ".pdf")) %>% ggsave(plot = p)
}

# All in one
pml_omics_tab <- tibble::tribble(
  ~omics_type, ~norm_method, ~reduction_model, ~fpath,
  "pml_citeseq", "SCT", "wnn.umap", "CITE_seq/objects/integrated/pbmc.cite_seq.integrated.pca_umap_clustered.annotated.rds",
  # "pml_multiome", "SCT", "umap", "Multiome/objects/integrated/pbmc.multiome.integrated.pca_umap_clustered.annotated.rds",
  "pml_rnaseq", "SCT", "umap", "RNA_seq/objects/integrated/pbmc.rna_seq.integrated.pca_umap_clustered.annotated.rds",
  # "pml_atacseq", "SCT", "umap", "ATAC_seq/objects/integrated/pbmc.atac_seq.integrated.pca_umap_clustered.annotated.rds",
)

selected_donors <- c(
  "PML0002", "PML0009", "PML0017", "PML0022", "PML0025", # CITE-seq
  "PML0033", "PML0055", "PML0058", "PML0060", "PML0063", "PML0066", "PML0070" # RNA-seq
)
# Map selected sc-omics results to reference.
pbmc_list <- apply(pml_omics_tab, 1, function(vec, .pbmc_ref, .overwrite) {
  # vec <- c(omics_type = "pml_rnaseq", norm_method = "SCT", reduction_model = "umap", fpath = "RNA_seq/objects/integrated/pbmc.rna_seq.integrated.pca_umap_clustered.annotated.rds")
  fpath <- vec["fpath"]
  red_model <- vec["reduction_model"]
  omics_type <- vec["omics_type"]
  norm_method <- vec["norm_method"]

  save_to <- file.path(object_dir, paste0(omics_type, ".rds"))
  if (!file.exists(save_to) || .overwrite) {
    cat("[I]: Mapping", omics_type, "\n")
    per_omics <- file.path(proj_dir, "outputs/analysis", fpath) %>% readRDS()
    per_omics <- per_omics[, per_omics$Vireo_assignment %in% selected_donors]
    anchors <- FindTransferAnchors(.pbmc_ref, per_omics, normalization.method = norm_method, reference.reduction = "spca", dims = 1:50, verbose = FALSE)
    per_omics <- MapQuery(anchors, per_omics, .pbmc_ref, refdata = list(celltype.l1 = "celltype.l1", celltype.l2 = "celltype.l2"), reference.reduction = "spca", verbose = FALSE)
    per_omics <- RunUMAP(per_omics, reduction = "ref.spca", dims = 1:30, verbose = FALSE)
    per_omics <- FindNeighbors(per_omics, reduction = "ref.spca", verbose = FALSE)
    per_omics <- FindClusters(per_omics, verbose = FALSE)
    cat("[I]: Dumping into disk...", omics_type, "\n")
    saveRDS(per_omics, save_to)
  } else {
    cat("[I]: Loading from disk...", omics_type, "\n")
    per_omics <- readRDS(save_to)
  }

  for (per_idents in c(paste0("predicted.celltype.", c("l1", "l2")))) {
    plot_width <- c("predicted.celltype.l1" = 6, "predicted.celltype.l2" = 10)[per_idents]
    plot_height <- c("predicted.celltype.l1" = 5, "predicted.celltype.l2" = 7)[per_idents]

    p <- DimPlot(per_omics, reduction = "ref.umap", group.by = per_idents, label = TRUE, repel = TRUE, pt.size = 0.75)
    save_to <- file.path(plot_dir, paste0(omics_type, ".ref_umap.", per_idents, ".pdf"))
    ggsave(save_to, plot = p, width = plot_width, height = plot_height)

    p <- DimPlot(per_omics, reduction = "umap", group.by = per_idents, label = TRUE, repel = TRUE, pt.size = 0.75)
    save_to <- file.path(plot_dir, paste0(omics_type, ".umap.", per_idents, ".pdf"))
    ggsave(save_to, plot = p, width = plot_width, height = plot_height)
  }

  return(invisible(per_omics))
}, .overwrite = TRUE, .pbmc_ref = pbmc_reference) %>%
  purrr::set_names(pml_omics_tab$omics_type)
