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

chrom_xy_genes <- fread(file.path(proj_dir, "inputs/references/sextual_chromosome_genes.txt"), header = F) %>% pull(V2)


# Load GEX reference dataset from Hao etal 2020
pbmc_reference <- file.path(proj_dir, "inputs/references/Hao_etal_2020/pbmc_multimodal.h5seurat") %>% LoadH5Seurat()
Misc(pbmc_reference[["wnn.umap"]], slot = "model")$num_precomputed_nns <- 1
for (per_idents in c(paste0("celltype.", c("l1", "l2")))) {
  save_to <- file.path(plot_dir, paste0("pbmc_reference_rnaseq.umap.", per_idents, ".pdf"))
  if (!file.exists(save_to)) {
    Idents(pbmc_reference) <- per_idents
    p <- DimPlot(pbmc_reference, reduction = "wnn.umap", label = TRUE, repel = TRUE) + NoLegend()
    ggsave(save_to, plot = p)
  } else {
    cat("[I]: Skip", save_to, "\n")
  }
}

# All in one
pml_omics_tab <- tibble::tribble(
  ~omics_type, ~norm_method, ~reduction_model, ~fpath,
  "pml_citeseq", "SCT", "wnn.umap", "CITE_seq/objects/integrated/pbmc.cite_seq.integrated.pca_umap_clustered.annotated.rds",
  "pml_rnaseq", "SCT", "umap", "RNA_seq/objects/integrated/pbmc.rna_seq.integrated.pca_umap_clustered.annotated.rds",
  # "pml_multiome", "SCT", "umap", "Multiome/objects/integrated/pbmc.multiome.integrated.pca_umap_clustered.annotated.rds",
  # "pml_atacseq", "SCT", "umap", "ATAC_seq/objects/integrated/pbmc.atac_seq.integrated.pca_umap_clustered.annotated.rds",
)

selected_donors <- c(
  "PML0002", "PML0009", "PML0017", "PML0022", "PML0025",
  "PML0013", "PML0033", "PML0055", "PML0058", "PML0060", "PML0061", "PML0063", "PML0070"
)
# Map selected sc-omics results to reference.
pbmc_list <- apply(pml_omics_tab, 1, function(vec, .pbmc_ref, .overwrite) {
  fpath <- vec["fpath"]
  omics_type <- vec["omics_type"]
  reduc_model <- vec["reduction_model"]
  norm_method <- vec["norm_method"]

  save_to <- file.path(object_dir, paste0(omics_type, ".rds"))
  if (!file.exists(save_to) || .overwrite) {
    cat("[I]: Mapping", omics_type, "\n")
    per_omics <- file.path(proj_dir, "outputs/analysis", fpath) %>% readRDS()
    per_omics <- per_omics[, per_omics$Vireo_assignment %in% selected_donors]
    anchors <- FindTransferAnchors(.pbmc_ref, per_omics, normalization.method = norm_method, reference.reduction = "spca", query.assay = "integrated", dims = 1:50, verbose = FALSE)
    per_omics <- MapQuery(anchors, per_omics, .pbmc_ref, refdata = list(celltype.l1 = "celltype.l1", celltype.l2 = "celltype.l2"), reference.reduction = "spca", verbose = FALSE)
    per_omics <- RunUMAP(per_omics, reduction = "ref.spca", reduction.name = "ref.reumap", dims = 1:30, verbose = FALSE)
    per_omics <- FindNeighbors(per_omics, reduction = "ref.spca", verbose = FALSE)
    per_omics <- FindClusters(per_omics, verbose = FALSE)
    per_omics@meta.data <- per_omics@meta.data %>%
      as.data.frame() %>%
      dplyr::mutate(predicted.celltype.l2 = dplyr::case_when(
        predicted.celltype.l2 == "CD8 Proliferating" ~ "CD8 Naive",
        predicted.celltype.l2 == "CD8 Proliferating" ~ "CD4 Navie",
        predicted.celltype.l2 %in% c("NK Proliferating", "NK_CD56bright") ~ "NK",
        predicted.celltype.l2 %in% c("ASDC", "MAIT", "ILC", "HSPC", "Eryth", "Plasmablast", "Platelet") ~ "other",
        TRUE ~ predicted.celltype.l2
      )) %>%
      dplyr::mutate(predicted.celltype.l1 = dplyr::case_when(
        predicted.celltype.l2 %in% c("CD4 CTL", "CD4 Navie", "CD4 Proliferating", "CD4 TCM", "CD4 TEM", "Treg") ~ "CD4 T",
        predicted.celltype.l2 %in% c("CD8 TEM", "CD8 Naive", "CD8 Proliferating", "CD8 TCM") ~ "CD8 T",
        predicted.celltype.l2 %in% c("CD14 Mono", "CD16 Mono") ~ "Mono",
        predicted.celltype.l2 %in% c("B intermediate", "B memory", "B naive") ~ "B",
        predicted.celltype.l2 %in% c("NK", "NK Proliferating", "NK_CD56bright") ~ "NK",
        predicted.celltype.l2 %in% c("cDC2", "pDC") ~ "DC",
        predicted.celltype.l2 %in% c("dnT", "gdT") ~ "other T",
        TRUE ~ "other"
      ))

    cat("[I]: Dumping into disk...", omics_type, "\n")
    saveRDS(per_omics, save_to)
  } else {
    cat("[I]: Loading from disk...", omics_type, "\n")
    per_omics <- readRDS(save_to)
  }

  for (per_idents in c("seurat_clusters", paste0("predicted.celltype.", c("l1", "l2")), "Vireo_assignment", "Response", "Timepoint", "Sex")) {
    plot_width <- c("predicted.celltype.l1" = 6, "predicted.celltype.l2" = 10)[per_idents]
    plot_height <- c("predicted.celltype.l1" = 5, "predicted.celltype.l2" = 7)[per_idents]

    save_to <- file.path(plot_dir, paste0(omics_type, ".ref_umap.", per_idents, ".pdf"))
    p <- DimPlot(per_omics, reduction = "ref.reumap", group.by = per_idents, label = TRUE, repel = TRUE, pt.size = 0.1)
    ggsave(save_to, plot = p, width = plot_width, height = plot_height)

    save_to <- file.path(plot_dir, paste0(omics_type, ".umap.", per_idents, ".pdf"))
    p <- DimPlot(per_omics, reduction = reduc_model, group.by = per_idents, label = TRUE, repel = TRUE, pt.size = 0.1)
    ggsave(save_to, plot = p, width = plot_width, height = plot_height)
  }

  return(invisible(per_omics))
}, .overwrite = FALSE, .pbmc_ref = pbmc_reference) %>%
  purrr::set_names(pml_omics_tab$omics_type)
