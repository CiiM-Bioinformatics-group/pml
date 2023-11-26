#!/usr/bin/env Rscript
# File: cite_seq.r
# Author: Zhenhua Zhang
# E-mail: zhenhua.zhang217@gmail.com
# Created: Jun 06, 2023
# Updated: Nov 22, 2023

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
})


projdir <- "~/Documents/projects/wp_pml"
deg_dir <- file.path(projdir, "outputs/analysis/CITE_seq/deg")
plot_dir <- file.path(projdir, "outputs/analysis/CITE_seq/plots")
object_dir <- file.path(projdir, "outputs/analysis/CITE_seq/objects")


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
meta_data_path <- file.path(projdir, "misc/metadata/clinical_parameters.csv")
major_meta_tab <- fread(meta_data_path) %>%
  dplyr::select(Patient_ID, Age, Sex, BMI, Response) %>%
  dplyr::mutate(Patient_ID = stringr::str_remove(Patient_ID, "Cure-"))

# Sequencing meta data, samples per pool
other_meta_tab <- tibble::tribble(
  ~Patient_ID, ~pool, ~Timepoint, # ~Timepoint_raw,
  "PML0002", "CITEpool2", "3M", # "3 months",
  "PML0002", "CITEpool4", "BL", # "BL",
  "PML0009", "CITEpool3", "3M", # "3 months",
  "PML0009", "CITEpool2", "BL", # "BL",
  "PML0017", "CITEpool1", "6W", # "6 weeks",
  "PML0017", "CITEpool4", "BL", # "BL",
  "PML0022", "CITEpool1", "3M", # "3 months",
  "PML0022", "CITEpool2", "6W", # "6 weeks",
  "PML0022", "CITEpool3", "BL", # "BL",
  "PML0025", "CITEpool1", "6M", # "6 months",
  "PML0025", "CITEpool4", "BL", # "BL",
)


# Load data and basic QC
selected_pools <- 1:4
pools <- paste0("CITEpool", selected_pools) %>% purrr::set_names(.)
bcmat_path <- file.path(projdir, "outputs/readcounts/CITE_seq", pools, "outs") %>% purrr::set_names(pools)
dmres_path <- file.path(projdir, "outputs/demultiplex/CITE_seq", pools, "vireo_outdir_ref/donor_ids.tsv") %>% purrr::set_names(pools)

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
      dplyr::mutate(pool = pn) %>%
      dplyr::left_join(dm_tab, by = "Cellbarcode") %>%
      dplyr::left_join(major_meta_tab, by = c("Vireo_assignment" = "Patient_ID")) %>%
      dplyr::left_join(other_meta_tab, by = c("Vireo_assignment" = "Patient_ID", "pool"))
    po@meta.data <- po@meta.data %>%
      dplyr::mutate(Cellbarcode = rownames(.) %>% stringr::str_remove("_[0-9]$"), pool = pn) %>%
      dplyr::left_join(extra_meta_tab, by = c("pool", "Cellbarcode")) %>%
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
}, .dmres_path = dmres_path, .bcmat_path = bcmat_path, .overwrite = TRUE)


# Filtering cells using soft threshold quantile(probs = c(0.05, 0.95)) on each parameter.
soft_threshold <- tibble::tribble(
  ~pool, ~min_ncount_rna, ~max_ncount_rna, ~min_nfeature_rna, ~max_nfeature_rna, ~max_percent_mt,
  "CITEpool1", 2087, 17112, 986, 4190, 9.4,
  "CITEpool2", 1811, 11152, 773, 3380, 9.8,
  "CITEpool3", 1763, 20000, 574, 3742, 9.5,
  "CITEpool4", 2037, 18513, 1041, 4393, 9.4,
)

selected_donors <- c("PML0002", "PML0009", "PML0017", "PML0022", "PML0025")
pbmc <- lapply(pools, function(pn, .pbmc, .hard_threshold, .cc_genes) {
  params <- .hard_threshold %>% dplyr::filter(pool == pn)
  min_ncount_rna <- params["min_ncount_rna"] %>% as.integer()
  max_ncount_rna <- params["max_ncount_rna"] %>% as.integer()
  min_nfeature_rna <- params["min_nfeature_rna"] %>% as.integer()
  max_nfeature_rna <- params["max_nfeature_rna"] %>% as.integer()
  max_percent_mt <- params["max_percent_mt"] %>% as.double()

  po <- .pbmc[[pn]]

  po$pool <- pn
  tar_cells <- po@meta.data %>%
    as.data.frame() %>%
    dplyr::filter(min_ncount_rna <= nCount_RNA, nCount_RNA <= max_ncount_rna, min_nfeature_rna <= nFeature_RNA, nFeature_RNA <= max_nfeature_rna, 100 <= nFeature_ADT, 500 <= nCount_ADT, nCount_ADT <= 20000, percent_mt < max_percent_mt, percent_rb < 50) %>%
    dplyr::filter(Vireo_assignment %in% selected_donors) %>%
    dplyr::filter(!(Vireo_assignment == "PML0002" & Timepoint == "3M")) %>%
    rownames()

  tar_features <- c(rownames(po@assays$RNA), rownames(po@assays$ADT)) %>% purrr::discard(~stringr::str_detect(.x, "^MT-|^RP[LS]"))

  po <- po[tar_features, tar_cells] %>%
    NormalizeData(assay = "RNA") %>%
    CellCycleScoring(assay = "RNA", s.features = .cc_genes$s.genes, g2m.features = .cc_genes$g2m.genes, set.ident = TRUE) %>%
    SCTransform(assay = "RNA", vst.flavor = "v2", vars.to.regress = c("percent_mt", "S.Score", "G2M.Score"), method = "glmGamPoi", variable.features.n = 3000) %>%
    NormalizeData(assay = "ADT", normalization.method = "CLR", margin = 2) %>%
    RunPCA(assay = "SCT", npcs = 30, verbose = FALSE)

  invisible(po)
}, .pbmc = pbmc_list, .hard_threshold = soft_threshold, .cc_genes = cc.genes)

# Integrate multiple pools
features <- SelectIntegrationFeatures(pbmc)
pbmc <- PrepSCTIntegration(pbmc, anchor.features = features)
anchors <- FindIntegrationAnchors(pbmc, normalization.method = "SCT", anchor.features = features, reduction = "rpca", dims = 1:30, k.anchor = 20)
pbmc_int <- IntegrateData(anchors, normalization.method = "SCT", dims = 1:30)

# Do PCA and regress out cell cycle scores
pbmc_int <- RunPCA(pbmc_int, assay = "integrated", npcs = 30, reduction.name = "int_pca") %>%
  ScaleData(assay = "ADT", vars.to.regress = c("percent_mt", "S.Score", "G2M.Score")) %>%
  RunPCA(assay = "ADT", features = rownames(pbmc_int[["ADT"]]), npcs = 30, reduction.name = "adt_pca")

# Find neighbors using multiple modal data and do UMAP
pbmc_int <- FindMultiModalNeighbors(pbmc_int, reduction.list = list("int_pca", "adt_pca"), dims.list = list(1:30, 1:30), modality.weight.name = "RNA.weight")
pbmc_int <- RunUMAP(pbmc_int, nn.name = "weighted.nn", reduction = "int_pca", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_", dims = 1:30)
pbmc_int <- FindClusters(pbmc_int, graph.name = "wsnn", resolution = 0.1)

# Prepare for DEG analysis
pbmc_int <- PrepSCTFindMarkers(pbmc_int)

# Annotate the cells by Azimuth
DefaultAssay(pbmc_int) <- "integrated"
pbmc_int <- RunAzimuth(pbmc_int, reference = "pbmcref")

# Save the data to disk
save_to <- file.path(object_dir, "integrated/pbmc.cite_seq.integrated.pca_umap_clustered.annotated.rds")
saveRDS(pbmc_int, save_to)


#
## Result overview
#
# Estimate the entropy by CellMixS
pbmc_cmb_sce <- as.SingleCellExperiment(pbmc_int)
reducedDim(pbmc_cmb_sce, "UMAP") <- pbmc_int@reductions$wnn.umap@cell.embeddings
reducedDim(pbmc_cmb_sce, "PCA") <- pbmc_int@reductions$int_pca@cell.embeddings
pbmc_cmb_sce <- entropy(pbmc_cmb_sce, "pool", k = 20)
save_to <- file.path(plot_dir, "integrated/pbmc.cite_seq.integrated.umap_entropy.pdf")
pdf(save_to, width = 12, height = 5)
visOverview(pbmc_cmb_sce, "pool", dim_red = "UMAP", metric = "entropy")
dev.off()

# Check the clusters, PCA
for (ppca in c("int_pca", "adt_pca")) {
  p <- DimPlot(pbmc_int, reduction = ppca)
  save_to <- file.path(plot_dir, "integrated", paste0("pbmc.cite_seq.integrated.", ppca, ".pdf"))
  ggsave(save_to, plot = p)
}

# Check the clusters, UMAP
for (per_group in c("seurat_clusters", "predicted.celltype.l1", "predicted.celltype.l2", "pool", "Response", "Sex", "Vireo_assignment", "Timepoint")) {
  p <- DimPlot(pbmc_int, reduction = "wnn.umap", group.by = per_group)
  nm <- stringr::str_replace_all(per_group, "\\.", "_")
  save_to <- file.path(plot_dir, "integrated", paste0("pbmc.cite_seq.integrated.umap_by_", nm, ".pdf"))
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

# Check example
p1 <- FeaturePlot(pbmc_int, "Hu.CD16", cols = c("lightgrey", "darkgreen"), reduction = "wnn.umap") + ggtitle("CD16 protein")
p2 <- FeaturePlot(pbmc_int, "FCGR3A", reduction = "wnn.umap") + ggtitle("FCGR3A(CD16) expression")
p <- p1 | p2
save_to <- file.path(plot_dir, "integrated/pbmc.cite_seq.integrated.CD16_gex.feature_plot.pdf")
ggsave(save_to, plot = p, width = 12, height = 6)


#
## Cell proportion
#
cell_propotion_tab <- pbmc_int@meta.data %>%
  dplyr::group_by(Vireo_assignment, Timepoint, Response, predicted.celltype.l1) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::mutate(prop = n / sum(n)) %>%
  dplyr::arrange(desc(prop)) %>%
  dplyr::mutate(predicted.celltype.l1 = factor(predicted.celltype.l1, levels = c("Mono", "CD8 T", "CD4 T", "B", "NK", "DC", "other T", "other"))) %>%
  dplyr::mutate(Timepoint = factor(Timepoint, levels = c("BL", "6W", "3M"))) %>%
  dplyr::ungroup() %>%
  dplyr::select(PatientID = Vireo_assignment, Timepoint, Response, predicted.celltype.l1, prop)
save_to <- file.path(plot_dir, "cell_proportion/pbmc.cite_seq.integrated.cell_propotion.csv")
fwrite(cell_propotion_tab, save_to)

# 1. Cell proportion overview
p <- cell_propotion_tab %>%
  ggplot() +
  geom_line(aes(x = Timepoint, y = prop, group = PatientID)) +
  geom_point(aes(x = Timepoint, y = prop, color = Response), alpha = 0.75) +
  scale_fill_npg() +
  theme_classic() +
  theme(axis.text.y = element_text(size = 8), axis.title = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(y = "Proportion", x = "Cell type") +
  facet_grid(~predicted.celltype.l1, space = "free_x")
save_to <- file.path(plot_dir, "cell_proportion/pbmc.cite_seq.integrated.cell_proportion.pdf")
ggsave(save_to, plot = p, width = 8.5, height = 3)

# 2. Cell proportion by patient
p <- cell_propotion_tab %>%
  ggplot() +
  geom_point(aes(x = PatientID, y = prop, color = predicted.celltype.l1)) +
  scale_color_npg() +
  theme_classic() +
  theme(axis.text.y = element_text(size = 8), axis.title = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(y = "Proportion", x = "Cell type") +
  facet_grid(~Timepoint, scales = "free_x", space = "free_x")
save_to <- file.path(plot_dir, "cell_proportion/pbmc.cite_seq.integrated.cell_proportion_by_patient.pdf")
ggsave(save_to, plot = p, width = 8.5, height = 3)

# 3. Cell proportion by timepoints
p <- cell_propotion_tab %>%
  ggplot() +
  geom_col(aes(x = Timepoint, y = prop, fill = predicted.celltype.l1), width = 0.985) +
  scale_fill_npg() +
  theme_classic() +
  theme(axis.text.y = element_text(size = 8), axis.title = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(y = "Proportion", x = "Cell type") +
  facet_grid(~PatientID, scales = "free_x", space = "free_x")
save_to <- file.path(plot_dir, "cell_proportion/pbmc.cite_seq.integrated.cell_proportion_by_timepoints.pdf")
ggsave(save_to, plot = p, width = 8.5, height = 3)


#
## DE gene analysis, using unnormalized RNA data
#
de_tab <- NULL
# 1. DE genes between responder and non-responder per cell type per time point
Idents(pbmc_int) <- paste(pbmc_int$predicted.celltype.l1, pbmc_int$Response, pbmc_int$Timepoint, sep = "_")
for (pct in tar_cell_types) {
  for (ptp in c("BL", "6W", "3M")) {
    id_1 <- paste0(pct, "_Responder_", ptp)
    id_2 <- paste0(pct, "_Non-responder_", ptp)
    tryCatch({
      de_tab <- FindMarkers(pbmc_int, ident.1 = id_1, ident.2 = id_2, min.pct = 0.2, logfc.threshold = 0.2) %>% 
        dplyr::mutate(gene_symbol = rownames(.), donors = "All", celltype = pct, comparison = paste0("Rs.vs.NRs_", ptp)) %>%
        dplyr::mutate(direction = dplyr::if_else(avg_log2FC > 0, "up", "dw")) %>%
        rbind(de_tab, .)
    }, error = function(e) cat(e$message, "; ", paste(id_1, id_2, "failed, comparison 1\n")))
  }
}

# 2. DE genes between three time points in responders per cell type.
Idents(pbmc_int) <- paste(pbmc_int$predicted.celltype.l1, pbmc_int$Timepoint, sep = "_")
tar_cells <- pbmc_int@meta.data %>% dplyr::filter(Response == "Non-responder") %>% rownames()
for (pct in tar_cell_types) {
  id_1 <- paste0(pct, "_BL")
  id_2 <- paste0(pct, "_6W")
  id_3 <- paste0(pct, "_3M")
  tryCatch({
    de_tab <- FindMarkers(pbmc_int[, tar_cells], ident.1 = id_3, ident.2 = id_1, min.pct = 0.2, logfc.threshold = 0.2) %>%
      dplyr::mutate(gene_symbol = rownames(.), donors = "Responder", celltype = pct, comparison = "3M.vs.BL") %>%
      dplyr::mutate(direction = dplyr::if_else(avg_log2FC > 0, "up", "dw")) %>%
      rbind(de_tab)
  }, error = function(e) cat(e$message, "; ", paste(id_2, id_1, "failed, comparison 2\n")))
  tryCatch({
    de_tab <- FindMarkers(pbmc_int[, tar_cells], ident.1 = id_2, ident.2 = id_1, min.pct = 0.2, logfc.threshold = 0.2) %>%
      dplyr::mutate(gene_symbol = rownames(.), donors = "Responder", celltype = pct, comparison = "6W.vs.BL") %>%
      dplyr::mutate(direction = dplyr::if_else(avg_log2FC > 0, "up", "dw")) %>%
      rbind(de_tab, .)
  }, error = function(e) cat(e$message, "; ", paste(id_2, id_1, "failed, comparison 2\n")))
  tryCatch({
    de_tab <- FindMarkers(pbmc_int[, tar_cells], ident.1 = id_3, ident.2 = id_2, min.pct = 0.2, logfc.threshold = 0.2) %>%
      dplyr::mutate(gene_symbol = rownames(.), donors = "Responder", celltype = pct, comparison = "3M.vs.6w") %>%
      dplyr::mutate(direction = dplyr::if_else(avg_log2FC > 0, "up", "dw")) %>%
      rbind(de_tab, .)
  }, error = function(e) cat(e$message, "; ", paste(id_2, id_1, "failed, comparison 2\n")))
}

# 3. DE genes between three time points in non-responders per cell type.
Idents(pbmc_int) <- paste(pbmc_int$predicted.celltype.l1, pbmc_int$Timepoint, sep = "_")
tar_cells <- pbmc_int@meta.data %>% dplyr::filter(Response == "Non-responder") %>% rownames()
for (pct in tar_cell_types) {
  id_1 <- paste0(pct, "_BL")
  id_2 <- paste0(pct, "_3M")
  tryCatch({
    de_tab <- FindMarkers(pbmc_int[, tar_cells], ident.1 = id_2, ident.2 = id_1, min.pct = 0.2, logfc.threshold = 0.2) %>%
      dplyr::mutate(gene_symbol = rownames(.), donors = "Non-responder", celltype = pct, comparison = "3m.vs.BL") %>%
      dplyr::mutate(direction = dplyr::if_else(avg_log2FC > 0, "up", "dw")) %>%
      rbind(de_tab, .)
  }, error = function(e) cat(e$message, "; ", paste(id_2, id_1, "failed, comparison 3\n")) )
}

save_to <- file.path(projdir, "outputs/analysis/CITE_seq/deg/pbmc.cite_seq.integrated.de_gene.csv")
fwrite(de_tab, save_to)
