#!/usr/bin/env Rscript
# File: multiome.r
# Author: Zhenhua Zhang
# E-mail: zhenhua.zhang217@gmail.com
# Created: Jun 12, 2023
# Updated: Nov 23, 2023

suppressPackageStartupMessages({
  library(lobstr)
  library(tidyverse)
  library(data.table)
  library(ggsci)

  library(Seurat)
  library(Signac)
  library(Azimuth)
  library(CellMixS)
  library(SeuratDisk)
  library(EnsDb.Hsapiens.v86)
  library(SingleCellExperiment)
  library(BSgenome.Hsapiens.UCSC.hg38)
})

projdir <- "~/Documents/projects/wp_pml"
deg_dir <- file.path(projdir, "outputs/analysis/Multiome/deg")
plot_dir <- file.path(projdir, "outputs/analysis/Multiome/plots")
object_dir <- file.path(projdir, "outputs/analysis/Multiome/objects")


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
## Multiome analysis
#
meta_data_path <- file.path(projdir, "misc/metadata/clinical_parameters.csv")
major_meta_tab <- fread(meta_data_path) %>%
  dplyr::select(Patient_ID, Age, Sex, BMI, Response) %>%
  dplyr::mutate(Patient_ID = stringr::str_remove(Patient_ID, "Cure-"))

# Sequencing meta data, samples per pool
other_meta_tab <- tibble::tribble(
  ~Patient_ID, ~pool, ~Timepoint, # ~Timepoint_raw,
  "PML0002", "MOpool2", "3M", # "3 months",
  "PML0002", "MOpool4", "BL", # "BL",
  "PML0009", "MOpool3", "3M", # "3 months",
  "PML0009", "MOpool2", "BL", # "BL",
  "PML0017", "MOpool1", "6W", # "6 weeks",
  "PML0017", "MOpool4", "BL", # "BL",
  "PML0022", "MOpool1", "3M", # "3 months",
  "PML0022", "MOpool2", "6W", # "6 weeks",
  "PML0022", "MOpool3", "BL", # "BL",
  "PML0025", "MOpool1", "6M", # "6 months",
  "PML0025", "MOpool4", "BL", # "BL",
)

# Annotations
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"

# Load data and basic QC
selected_pools <- c(1, 3, 4)
pools <- paste0("MOpool", selected_pools) %>% purrr::set_names(.)
bcmat_path <- file.path(projdir, "outputs/readcounts/Multiome", pools, "outs") %>% purrr::set_names(pools)
dmres_path <- file.path(projdir, "outputs/demultiplex/Multiome", pools, "vireo_outdir_ref/donor_ids.tsv") %>% purrr::set_names(pools)

# Common features across pools
all_chroms <- paste0("chr", c(1:22, "X", "Y"))
common_features <- lapply(pools, function(pn, .bcmat_path, .overwrite = FALSE) {
  pp <- file.path(.bcmat_path[[pn]])
  peak_gr <- file.path(pp, "atac_peaks.bed") %>%
    read.table(col.names = c("chr", "start", "end")) %>%
    dplyr::filter(chr %in% all_chroms) %>%
    makeGRangesFromDataFrame()
}, .bcmat_path = bcmat_path, .overwrite = overwrite) %>%
  Reduce(c, .) %>%
  IRanges::reduce(ignore.strand = TRUE) %>%
  (function(gr) { w <- width(gr); gr[w < 2000 & w > 100]})

# Load data including GEX and CA
pbmc_list <- lapply(pools, function(pn, .overwrite, .dmres_path, .bcmat_path, .common_features, .annotations) { # pn <- pools[1]
  pp <- file.path(.bcmat_path[[pn]], "filtered_feature_bc_matrix")
  obj_save_to <- file.path(projdir, "outputs/analysis/Multiome/objects", paste0("pbmc.", pn, ".raw.rds"))
  if (!file.exists(obj_save_to) || .overwrite) {
    cat("[I]: Creating Seurat object for", pn, "\n")
    # Loading demultiplex results
    dm_tab <- fread(.dmres_path[pn]) %>%
      dplyr::select(Cellbarcode = cell, Vireo_assignment = donor_id, Vireo_prob_max = prob_max, Vireo_prob_doublet = prob_doublet, Vireo_dbl_llr = doublet_logLikRatio) %>%
      dplyr::mutate(Vireo_assignment = dplyr::case_when(Vireo_assignment %in% c("doublet", "unassigned") ~ Vireo_assignment, TRUE ~ paste0("PML", stringr::str_remove(Vireo_assignment, "PML") %>% stringr::str_pad(4, pad = "0"))))

    # Load CellRanger counts results
    pm <- Read10X(pp)

    # Create RNA assay
    po <- CreateSeuratObject(counts = pm$`Gene Expression`, project = "PMLMultiome")

    # MT and RB gene ratios
    po[["percent_mt"]] <- PercentageFeatureSet(po, pattern = "^MT-")
    po[["percent_rb"]] <- PercentageFeatureSet(po, pattern = "^RP[LS]")

    # Create ATAC assay
    atac_cells <- file.path(.bcmat_path[pn], "atac_fragments.tsv.gz") %>% CountFragments() %>% dplyr::filter(CB %in% colnames(po)) %>% dplyr::pull(CB)
    fragments <- CreateFragmentObject(path = file.path(.bcmat_path[pn], "atac_fragments.tsv.gz"), cells = atac_cells)
    counts <- FeatureMatrix(fragments = fragments, features = .common_features, cells = atac_cells)
    po[["ATAC"]] <- CreateChromatinAssay(counts = counts, sep = c(":", "-"), genome = "hg38", fragments = fragments, annotation = .annotations)

    # Add ATAC QC metrics
    po <- NucleosomeSignal(po, assay = "ATAC") %>% TSSEnrichment(fast = FALSE, assay = "ATAC")

    # Update meta data
    raw_barcodes <- colnames(po)
    extra_meta_tab <- read.csv(file.path(pp, "barcodes.tsv.gz"), col.names = "Cellbarcode") %>%
      dplyr::mutate(pool = pn) %>%
      dplyr::left_join(dm_tab, by = "Cellbarcode") %>%
      dplyr::left_join(major_meta_tab, by = c("Vireo_assignment" = "Patient_ID")) %>%
      dplyr::left_join(other_meta_tab, by = c("Vireo_assignment" = "Patient_ID", "pool"))
    po@meta.data <- po@meta.data %>%
      dplyr::mutate(Cellbarcode = rownames(.) %>% stringr::str_remove("_[0-9]$"), pool = pn) %>%
      dplyr::left_join(extra_meta_tab, by = c("pool", "Cellbarcode")) %>%
      (function(tab) {rownames(tab) <- raw_barcodes; tab})

    # QC plots
    p <- VlnPlot(object = po, features = c("nCount_RNA", "nFeature_RNA", "nCount_ATAC", "nFeature_ATAC", "percent_mt", "percent_rb"), ncol = 5, pt.size = 0)
    save_to <- file.path(plot_dir, "quality_control", paste("control_metrics", pn, "violin.pdf", sep = "."))
    ggsave(save_to, p, width = 12, height = 8)

    p <- FeatureScatter(po, feature1 = "nCount_RNA", feature2 = "percent_mt") +
      FeatureScatter(po, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
      FeatureScatter(po, feature1 = "nCount_ATAC", feature2 = "nFeature_ATAC")
    save_to <- file.path(plot_dir, "quality_control", paste("control_metrics", pn, "feature_scatter.pdf", sep = "."))
    ggsave(save_to, p, width = 12, height = 7)

    # Dump into disk
    cat("[I]: Pool", pn, "has", ncol(po), "cells and", nrow(po), "genes.\n")
    cat("[I]: Dumping into disk ...", pn, "\n"); saveRDS(po, obj_save_to)
  } else {
    cat("[I]: Loading Seurat object for", pn, "\n"); po <- readRDS(obj_save_to)
  }

  invisible(po)
}, .overwrite = TRUE, .dmres_path = dmres_path, .bcmat_path = bcmat_path, .common_features = common_features, .annotations = annotations)


# quantile(n[Count|Feature]_[RNA|ATAC], probs = c(0.05, 0.95))
hard_threshold <- tibble::tribble(
  ~pool, ~min_ncount_rna, ~max_ncount_rna, ~min_ncount_atac, ~max_ncount_atac, ~max_percent_mt,
  "MOpool1", 538, 3862, 835, 2369, 20,
  "MOpool2", 500, 25000, 750, 5000, 20,
  "MOpool3", 371, 5330, 537, 9914, 20,
  "MOpool4", 540, 4404, 1090, 4900, 20,
)

selected_donors <- c("PML0002", "PML0009", "PML0017", "PML0022", "PML0025")
pbmc <- lapply(pools, function(pn, .pbmc, .hard_threshold, .cc_genes) {
  params <- .hard_threshold %>% dplyr::filter(pool == pn)
  min_ncount_rna <- params["min_ncount_rna"] %>% as.integer()
  max_ncount_rna <- params["max_ncount_rna"] %>% as.integer()
  min_ncount_atac <- params["min_ncount_atac"] %>% as.integer()
  max_ncount_atac <- params["max_ncount_atac"] %>% as.integer()
  max_percent_mt <- params["max_percent_mt"] %>% as.double()

  po <- .pbmc[[pn]]
  po$pool <- pn

  tar_cells <- po@meta.data %>%
    as.data.frame() %>%
    dplyr::filter(min_ncount_rna <= nCount_RNA, nCount_RNA <= max_ncount_rna, min_ncount_atac <= nCount_ATAC, nCount_ATAC <= max_ncount_atac, percent_mt < max_percent_mt) %>%
    dplyr::filter(Vireo_assignment %in% selected_donors) %>%
    rownames()

  tar_features <- rownames(po@assays$RNA) %>% purrr::discard(~stringr::str_detect(.x, "^MT-|^RP[LS]"))

  po <- po[, tar_cells] %>%
    NormalizeData(assay = "RNA", verbose = FALSE) %>%
    CellCycleScoring(assay = "RNA", s.features = .cc_genes$s.genes, g2m.features = .cc_genes$g2m.genes, set.ident = TRUE, verbose = FALSE) %>%
    SCTransform(assay = "RNA", vst.flavor = "v2", vars.to.regress = c("percent_mt", "S.Score", "G2M.Score"), method = "glmGamPoi", variable.features.n = 3000, verbose = FALSE) %>%
    RunPCA(assay = "SCT", npcs = 30, verbose = FALSE)

  invisible(po)
}, .pbmc = pbmc_list, .hard_threshold = hard_threshold, .cc_genes = cc.genes)

# Integrate multiple pools
features <- SelectIntegrationFeatures(pbmc)
pbmc <- PrepSCTIntegration(pbmc, anchor.features = features)
anchors <- FindIntegrationAnchors(pbmc, normalization.method = "SCT", anchor.features = features, verbose = FALSE)
pbmc_int <- IntegrateData(anchors, normalization.method = "SCT", verbose = FALSE)

# Do PCA
pbmc_int <- RunPCA(pbmc_int, npcs = 30)

# Do UMAP and find clusters
pbmc_int <- RunUMAP(pbmc_int, dims = 1:30)
pbmc_int <- FindNeighbors(pbmc_int, dims = 1:30)
pbmc_int <- FindClusters(pbmc_int, resolution = 0.1)

# Prepare for DEG analysis
pbmc_int <- PrepSCTFindMarkers(pbmc_int)

# Annotate the cells by Azimuth
DefaultAssay(pbmc_int) <- "integrated"
pbmc_int <- RunAzimuth(pbmc_int, reference = "pbmcref")
# pbmc_int$percent.mt <- NULL

# Save the data to disk
save_to <- file.path(object_dir, "integrated/pbmc.multiome.integrated.pca_umap_clustered.annotated.rds")
saveRDS(pbmc_int, save_to)


#
## Result overview
#
# Estimate the entropy by CellMixS
pbmc_cmb_sce <- as.SingleCellExperiment(pbmc_int)
reducedDim(pbmc_cmb_sce, "UMAP") <- pbmc_int@reductions$umap@cell.embeddings
reducedDim(pbmc_cmb_sce, "PCA") <- pbmc_int@reductions$pca@cell.embeddings
pbmc_cmb_sce <- entropy(pbmc_cmb_sce, "pool", k = 20)
save_to <- file.path(plot_dir, "integrated/pbmc.multiome.integrated.umap_entropy.pdf")
pdf(save_to, width = 12, height = 5)
visOverview(pbmc_cmb_sce, "pool", dim_red = "UMAP", metric = "entropy")
dev.off()

# Check the clusters, PCA
for (ppca in c("pca")) {
  p <- DimPlot(pbmc_int, reduction = ppca)
  save_to <- file.path(plot_dir, "integrated", paste0("pbmc.multiome.integrated.", ppca, ".pdf"))
  ggsave(save_to, plot = p)
}

# Check the clusters, UMAP
for (per_group in c("seurat_clusters", "predicted.celltype.l1", "predicted.celltype.l2", "pool", "Response", "Sex", "Vireo_assignment", "Timepoint")) {
  p <- DimPlot(pbmc_int, group.by = per_group)
  nm <- stringr::str_replace_all(per_group, "\\.", "_")
  save_to <- file.path(plot_dir, "integrated", paste0("pbmc.multiome.integrated.umap_by_", nm, ".pdf"))
  plot_wd <- ifelse(per_group == "predicted.celltype.l2", 10, 7)
  ggsave(save_to, plot = p, width = plot_wd)
}

# Cell type markers
for (per_group in c("seurat_clusters", "predicted.celltype.l1", "predicted.celltype.l2")) {
  p <- DotPlot(pbmc_int, features = celltype_markers, group.by = per_group, cluster.idents = TRUE) + coord_flip() + RotatedAxis()
  nm <- stringr::str_replace_all(per_group, "\\.", "_")
  save_to <- file.path(plot_dir, "integrated", paste0("pbmc.cite_seq.celltype_markers.dotplot_by_", nm, ".pdf"))
  plot_wd <- pbmc_int[[per_group]] %>% unique() %>% nrow() %>% `*`(0.4) %>% `+`(2)
  ggsave(save_to, p, width = plot_wd, height = 12)
}
