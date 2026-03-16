#!/usr/bin/env Rscript
# File: coanalysis_with_public_data.r
# Author: Zhenhua Zhang
# E-mail: zhenhua.zhang217@gmail.com
# Created: Nov 13, 2025
# Updated: 

# Co-analysis with published dataset

options(stringsAsFactors = FALSE, data.table.verbose = FALSE, error = rlang::last_trace, future.globals.maxSize = 1024 * 10 ^ 10)
suppressPackageStartupMessages({
  library(lobstr)
  library(data.table)
  library(tidyverse)
  library(patchwork)

  library(Seurat)
  library(Azimuth)
  library(SeuratData)
  library(SeuratDisk)

  library(SingleCellExperiment)
  library(scDblFinder)
  library(RColorBrewer)
  library(harmony)

  library(clusterProfiler)
  library(ggVennDiagram)
  library(ggsci)

  library(DropletUtils) # Save the Seurat object into 10X read counts
  library(monocle3)
})


projdir <- "~/Documents/projects/wp_pml"
plot_dir <- file.path(projdir, "outputs/analysis/public_data/plots")
object_dir <- file.path(projdir, "outputs/analysis/public_data/objects")


#
## Loading and preprocess public data.
#
rb_gradient <- RColorBrewer::brewer.pal(11, "RdBu")
chrom_xy_genes <- fread("~/Documents/projects/wp_pml/inputs/references/sextual_chromosome_genes.txt", header = F) %>% pull(V2)


incl_hiv_pml <- FALSE
saving_flag <- ifelse(incl_hiv_pml, "incl_hiv_pml", "excl_hiv_pml")
save_to <- file.path(projdir, "outputs/analysis/public_data/objects", paste0("Deffner_etal_2024.", saving_flag, ".rds"))
if (file.exists(save_to)) {
  csf_pub <- readRDS(save_to)
} else {
  # Lymphoma related PML
  csf_lmp <- "inputs/public_data/scRNAseq_data_PML_patients_Deffner_et_al_CellReports_Medicine/Lymphoma_related_PML" %>%
    file.path(projdir, .) %>%
    Read10X() %>%
    CreateSeuratObject(count = ., project = "Lymphoma_related_PML", min.cells = 3, min.features = 200) %>%
    AddMetaData("Lymphoma_related_PML", "SequencingPool")
  DefaultAssay(csf_lmp) <- "RNA"
  csf_lmp[["percent_mt"]] <- PercentageFeatureSet(csf_lmp, pattern = "^MT-")
  csf_lmp[["percent_rb"]] <- PercentageFeatureSet(csf_lmp, pattern = "^RP[LS]")

  p <- VlnPlot(object = csf_lmp, features = c("nCount_RNA", "nFeature_RNA", "percent_mt", "percent_rb"), ncol = 5, pt.size = 0)
  save_to <- file.path(plot_dir, "Deffner_etal_2024.quality_control_metrics.Lymphoma_related_PML.violin.pdf")
  ggsave(save_to, p, width = 12, height = 8)

  csf_lmp$nCount_RNA %>% quantile(probs = c(0.05, 0.95)) %>% as.integer()
  csf_lmp$nFeature_RNA %>% quantile(probs = c(0.05, 0.95)) %>% as.integer()
  csf_lmp <- subset(csf_lmp, subset = 2613 <= nCount_RNA & nCount_RNA <= 8684 & 1275 <= nFeature_RNA & nFeature_RNA <= 2867 & percent_mt <= 10)

  if (incl_hiv_pml) {
    pca_on_assay <- "integrated"
    cluster_on_graph <- "integrated.nn"
    # HIV related PML
    csf_hiv <- "inputs/public_data/scRNAseq_data_PML_patients_Deffner_et_al_CellReports_Medicine/HIV_related_PML" %>%
      file.path(projdir, .) %>%
      Read10X() %>%
      CreateSeuratObject(count = ., project = "HIV_related_PML", min.cells = 3, min.features = 200) %>%
      AddMetaData("HIV_related_PML", "SequencingPool")

    DefaultAssay(csf_hiv) <- "RNA"
    csf_hiv[["percent_mt"]] <- PercentageFeatureSet(csf_hiv, pattern = "^MT-")
    csf_hiv[["percent_rb"]] <- PercentageFeatureSet(csf_hiv, pattern = "^RP[LS]")

    p <- VlnPlot(object = csf_hiv, features = c("nCount_RNA", "nFeature_RNA", "percent_mt", "percent_rb"), ncol = 5, pt.size = 0)
    save_to <- file.path(plot_dir, "Deffner_etal_2024.quality_control_metrics.HIV_related_PML.violin.pdf")
    ggsave(save_to, p, width = 12, height = 8)

    csf_hiv$nCount_RNA %>% quantile(probs = c(0.05, 0.95)) %>% as.integer()
    csf_hiv$nFeature_RNA %>% quantile(probs = c(0.05, 0.95)) %>% as.integer()
    csf_hiv <- subset(csf_hiv, subset = 1327 <= nCount_RNA & nCount_RNA <= 7047 & 874 <= nFeature_RNA & nFeature_RNA <= 2873 & percent_mt <= 10)

    # Check batch effect
    csf_merge <- merge(csf_lmp, csf_hiv)
    p <- VlnPlot(object = csf_merge, features = c("nCount_RNA", "nFeature_RNA", "percent_mt", "percent_rb"), ncol = 5, pt.size = 0)
    save_to <- file.path(plot_dir, "Deffner_etal_2024.quality_control_metrics.violin.pdf")
    ggsave(save_to, p, width = 12, height = 8)
      
    csf_merge <- NormalizeData(csf_merge, verbose = FALSE) %>%
      FindVariableFeatures(verbose = FALSE) %>%
      ScaleData(verbose = FALSE) %>%
      RunPCA(verbose = FALSE) %>%
      RunUMAP(dims = 1:30, verbose = FALSE)

    p <- DimPlot(csf_merge, reduction = "umap", group.by = "orig.ident") + NoLegend()
    save_to <- file.path(plot_dir, "Deffner_etal_2024.quality_control_metrics.umap.pdf")
    ggsave(save_to, p, width = 7, height = 7)

    # Integration
    csf_merge_list <- csf_merge %>% SplitObject(split.by = "orig.ident") %>%
      lapply(function(po) {
        tar_features <- c(rownames(po@assays$RNA)) %>% purrr::discard(~stringr::str_detect(.x, "^MT-"))
        po[tar_features, ] %>%
          NormalizeData(assay = "RNA", verbose = FALSE) %>%
          SCTransform(assay = "RNA", vst.flavor = "v2", vars.to.regress = c("percent_mt"), method = "glmGamPoi", variable.features.n = 2000, verbose = FALSE) %>%
          RunPCA(assay = "SCT", verbose = FALSE)
      })
    features <- SelectIntegrationFeatures(csf_merge_list, verbose = FALSE)
    csf_merge_list <- PrepSCTIntegration(csf_merge_list, anchor.features = features, verbose = FALSE)
    csf_pub <- FindIntegrationAnchors(csf_merge_list, normalization.method = "SCT", anchor.features = features, reduction = "rpca", dims = 1:30, k.anchor = 20, verbose = FALSE)
    csf_pub <- IntegrateData(csf_pub, normalization.method = "SCT", dims = 1:30, verbose = FALSE)
  } else {
    pca_on_assay <- "SCT"
    cluster_on_graph <- "SCT_nn"
    csf_pub <- csf_lmp %>%
      NormalizeData(verbose = FALSE) %>%
      SCTransform(assay = "RNA", vst.flavor = "v2", vars.to.regress = c("percent_mt"), method = "glmGamPoi", variable.features.n = 2000, verbose = FALSE)
  }

  # Do PCA and regress out cell cycle scores
  csf_pub <- RunPCA(csf_pub, assay = pca_on_assay, verbose = FALSE)
  csf_pub <- RunUMAP(csf_pub, dims = 1:15, assay = pca_on_assay, n.neighbors = 50, verbose = FALSE)
  csf_pub <- csf_pub %>% (function(obj) { # Define doublets
    DefaultAssay(obj) <- "SCT"
    csf_cmb_sce <- as.SingleCellExperiment(obj)
    reducedDim(csf_cmb_sce, "UMAP") <- obj@reductions$umap@cell.embeddings
    reducedDim(csf_cmb_sce, "PCA") <- obj@reductions$pca@cell.embeddings
    csf_cmb_sce <- scDblFinder(csf_cmb_sce, clusters = TRUE)
    obj$scDblFinder_class <- csf_cmb_sce$scDblFinder.class
    obj$scDblFinder_score <- csf_cmb_sce$scDblFinder.score

    obj
  })

  csf_pub <- FindNeighbors(csf_pub, verbose = FALSE)
  csf_pub <- FindClusters(csf_pub, graph.name = cluster_on_graph, resolution = 0.1, verbose = FALSE)
  csf_pub <- PrepSCTFindMarkers(csf_pub, verbose = FALSE)

  # Check integration
  for (per_group in c("seurat_clusters", "orig.ident", "nCount_RNA", "nFeature_RNA", "percent_rb", "scDblFinder_score")) {
    nm <- stringr::str_replace_all(per_group, "\\.", "_")
    umap_save_to <- file.path(plot_dir, paste0("Deffner_etal_2024.", saving_flag, ".umap_by_", nm, ".pdf"))
    if (per_group %in% c("seurat_clusters", "orig.ident")) {
      p <- DimPlot(csf_pub, reduction = "umap", group.by = per_group, pt.size = 1.0)
    } else {
      p <- FeaturePlot(csf_pub, reduction = "umap", features = per_group, pt.size = 1.0)
    }
    ggsave(umap_save_to, plot = p, width = 7, height = 7)
  }

  # Identify cell types, find markers
  csf_pub_deg_tbl <- FindAllMarkers(csf_pub, assay = "SCT", verbose = FALSE)
  marker_features <- csf_pub_deg_tbl %>%
    dplyr::filter(p_val_adj < 0.05, avg_log2FC > 0.75) %>%
    dplyr::group_by(cluster) %>%
    dplyr::slice_min(p_val_adj, n = 10, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::pull(gene) %>%
    unique()

  p <- DotPlot(csf_pub, features = marker_features) +
    scale_color_gradient2(low = rb_gradient[11], mid = rb_gradient[6], high = rb_gradient[1]) +
    coord_flip()
  plot_save_to <- file.path(plot_dir, paste0("Deffner_etal_2024.", saving_flag, ".celltype_markers.dotplot.pdf"))
  ggsave(plot_save_to, plot = p, width = 6, height = 12)

  # Predict cell types
  pbmc_ref_path <- "/home/zzhang/Documents/projects/resources/Azimuth/PBMC"
  csf_pub <- RunAzimuth(csf_pub, reference = pbmc_ref_path, assay = "SCT")
  DefaultAssay(csf_pub) <- "SCT"

  csf_pub@meta.data <- csf_pub@meta.data %>%
    dplyr::mutate(
      predicted.celltype.l2 = dplyr::case_when(
        predicted.celltype.l2 == "CD8 Proliferating" ~ "CD8 Naive",
        predicted.celltype.l2 == "CD8 Proliferating" ~ "CD4 Navie",
        predicted.celltype.l2 %in% c("NK Proliferating", "NK_CD56bright") ~ "NK",
        predicted.celltype.l2 %in% c("ASDC", "MAIT", "ILC", "HSPC", "Eryth", "Plasmablast", "Platelet") ~ "other",
        TRUE ~ predicted.celltype.l2
      ),
      predicted.celltype.l1 = dplyr::case_when(
        predicted.celltype.l2 %in% c("CD4 CTL", "CD4 Navie", "CD4 Proliferating", "CD4 TCM", "CD4 TEM", "Treg") ~ "CD4 T",
        predicted.celltype.l2 %in% c("CD8 TEM", "CD8 Naive", "CD8 Proliferating", "CD8 TCM") ~ "CD8 T",
        predicted.celltype.l2 %in% c("CD14 Mono", "CD16 Mono") ~ "Mono",
        predicted.celltype.l2 %in% c("B intermediate", "B memory", "B naive") ~ "B",
        predicted.celltype.l2 %in% c("NK", "NK Proliferating", "NK_CD56bright") ~ "NK",
        predicted.celltype.l2 %in% c("cDC2", "pDC") ~ "DC",
        predicted.celltype.l2 %in% c("dnT", "gdT") ~ "other T",
        TRUE ~ "other"
      ))

  umap_save_to <- file.path(plot_dir, paste0("Deffner_etal_2024.", saving_flag, ".umap_by_celltype.pdf"))
  p <- DimPlot(csf_pub, reduction = "umap", group.by = "predicted.celltype.l1", pt.size = 0.1)
  ggsave(umap_save_to, plot = p, width = 7, height = 7)

  umap_save_to <- file.path(plot_dir, paste0("Deffner_etal_2024.IFNG.", saving_flag, ".umap_by_celltype.pdf"))
  p <- FeaturePlot(csf_pub, reduction = "umap", features = "IFNG", pt.size = 0.1)
  ggsave(umap_save_to, plot = p, width = 7, height = 7)

  Idents(csf_pub) <- "predicted.celltype.l1"
  p <- DotPlot(csf_pub, features = marker_features) +
    RotatedAxis() +
    scale_color_gradient2(low = rb_gradient[11], mid = rb_gradient[6], high = rb_gradient[1]) +
    coord_flip()
  plot_save_to <- file.path(plot_dir, paste0("Deffner_etal_2024.", saving_flag, ".celltype_markers.dotplot.by_celltypes.pdf"))
  ggsave(plot_save_to, plot = p, width = 6, height = 12)

  saveRDS(csf_pub, save_to)
}


# Co-analysis with pml-citeseq
pml_int_save_to <- file.path(plot_dir, paste0("pml_integrated.coanalysis_with_public_data.rds"))
if (file.exists(pml_int_save_to)) {
  pml_int <- readRDS(pml_int_save_to)
} else {
  pbmc_cite <- "outputs/analysis/integrated/objects/pml_citeseq.rds" %>% file.path(projdir, .) %>% readRDS()
  DefaultAssay(pbmc_cite) <- "RNA"
  pbmc_cite[["ADT"]] <- NULL
  pbmc_cite[["SCT"]] <- NULL
  pbmc_cite[["integrated"]] <- NULL
  pbmc_cite[["wnn.umap"]] <- NULL
  pbmc_cite[["ref.reumap"]] <- NULL

  DefaultAssay(csf_pub) <- "RNA"
  #csf_pub[["SCT"]] <- NULL
  #csf_pub[["umap"]] <- NULL
  #csf_pub[["refAssay"]] <- NULL
  #csf_pub[["ref.umap"]] <- NULL

  var_to_regress <- c("percent_mt")
  pml_int_list <- pbmc_cite %>% SplitObject("SequencingPool") %>% c(list(csf_pub = csf_pub)) %>%
    lapply(function(o) {
      DefaultAssay(o) <- "RNA"
      SCTransform(o, assay = "RNA", vst.flavor = "v2", vars.to.regress = var_to_regress, method = "glmGamPoi", verbose = FALSE) %>%
        RunPCA(assay = "SCT", verbose = FALSE)
    })
  integration_features <- SelectIntegrationFeatures(pml_int_list, nfeatures = 3000)
  pml_int_list <- PrepSCTIntegration(pml_int_list, anchor.features = integration_features, verbose = FALSE)
  pml_int <- FindIntegrationAnchors(pml_int_list, normalization.method = "SCT", anchor.features = integration_features, reduction = "rpca", dims = 1:30, k.anchor = 20, verbose = FALSE)
  pml_int <- IntegrateData(pml_int, normalization.method = "SCT", dims = 1:30, verbose = FALSE)
  pml_int <- RunPCA(pml_int, assay = "integrated", verbose = FALSE)
  pml_int <- RunUMAP(pml_int, assay = "integrated", dims = 1:30, verbose = FALSE)
  pml_int <- FindNeighbors(pml_int, assay = "integrated", dims = 1:30, verbose = FALSE)
  pml_int <- FindClusters(pml_int, graph.name = "integrated_nn", resolution = 0.1, verbose = FALSE)
  pml_int <- PrepSCTFindMarkers(pml_int)

  cols_to_incl <- c(
    "orig.ident", "nCount_RNA", "nFeature_RNA", "percent_mt", "SequencingPool", "Vireo_assignment",
    "Vireo_prob_max", "Vireo_prob_doublet", "Vireo_dbl_llr", "Age", "Sex", "BMI", "Response", "Timepoint", "nCount_SCT",
    "nFeature_SCT", "scDblFinder_class", "scDblFinder_score", "predicted.celltype.l1.score", "predicted.celltype.l1",
    "predicted.celltype.l2.score", "predicted.celltype.l2", "integrated_nn_res.0.1"
  )
  pml_int@meta.data <- pml_int@meta.data %>%
    dplyr::select(dplyr::all_of(cols_to_incl)) %>%
    dplyr::mutate(
      Timepoint = dplyr::if_else(is.na(Timepoint), "BL", Timepoint),
      Age = dplyr::if_else(is.na(Age), 70.5, Age),
      Sex = dplyr::if_else(is.na(Sex), "<NA>", Sex),
      BMI = dplyr::if_else(is.na(BMI), "<NA>", BMI),
      Response = dplyr::if_else(is.na(Response), "Unknown", Response),
      Vireo_assignment = dplyr::if_else(is.na(Vireo_assignment), "O-PML", Vireo_assignment),
      predicted.celltype.l1 = factor(predicted.celltype.l1, levels = c("CD8 T", "CD4 T", "NK", "Mono", "B", "DC", "other T", "other"))
    )

  # Check integration results.
  umap_save_to <- file.path(plot_dir, paste0("pml_integrated.umap_by_celltype.pdf"))
  p <- DimPlot(pml_int, reduction = "umap", group.by = "predicted.celltype.l1", pt.size = 0.1)
  ggsave(umap_save_to, plot = p, width = 7.25, height = 7)

  umap_save_to <- file.path(plot_dir, paste0("pml_integrated.umap_by_sequencingpool.pdf"))
  p <- DimPlot(pml_int, reduction = "umap", group.by = "SequencingPool", pt.size = 0.1)
  ggsave(umap_save_to, plot = p, width = 8.75, height = 7)

  umap_save_to <- file.path(plot_dir, paste0("pml_integrated.umap_by_datasource.pdf"))
  p <- DimPlot(pml_int, reduction = "umap", group.by = "orig.ident", pt.size = 0.1)
  ggsave(umap_save_to, plot = p, width = 8.75, height = 7)

  umap_save_to <- file.path(plot_dir, paste0("pml_integrated.umap_by_timepoint.pdf"))
  p <- DimPlot(pml_int, reduction = "umap", group.by = "Timepoint", pt.size = 0.1)
  ggsave(umap_save_to, plot = p, width = 7.25, height = 7)

  # Save to disk
  # pml_int[["prediction.score.celltype.l3"]] <- NULL
  # pml_int[["refAssay"]] <- NULL
  saveRDS(pml_int, pml_int_save_to)
}


# Check features
# CCL5, CXCL10 (inflammation chemokines)
# CCL2, CCL4
# CCR2, CCR5, CXCR3
# ITGA4, STXBP2, LY9
pbmc_cite <- "outputs/analysis/integrated/objects/pml_citeseq.rds" %>% file.path(projdir, .) %>% readRDS()
DefaultAssay(pbmc_cite) <- "RNA"
Idents(pbmc_cite) <- "predicted.celltype.l1"

Deffner_etal_2024_makers <- c(
  "CCL2", "CCL3", "CCL4", "CCL5", "CCL20",
  "CXCL1", "CXCL10", "CXCL12", "CXCL13",
  "CCR2", "CCR5", "CXCR3", "VCAM1",
  "ITGA4", "STXBP2", "LY9"
)

selected_cells <- pbmc_cite@meta.data %>% dplyr::filter(Timepoint %in% "BL") %>% rownames()
umap_save_to <- file.path(plot_dir, paste0("pbmc_cite.umap_by_Deffner_etal_2024_makers.pdf"))
p <- (FeaturePlot(pbmc_cite[, selected_cells], reduction = "wnn.umap", features = Deffner_etal_2024_makers, pt.size = 0.5) &
  scale_color_gradient(low = "gray95", high = rb_gradient[1], limits = c(0, 5.0)) &
  labs(x = NULL, y = NULL, color = "Average\nExpression") &
  theme(title = element_text(size = 11), axis.text = element_text(size = 9))) + plot_layout(guides = "collect") 
ggsave(umap_save_to, plot = p, width = 8.25, height = 7.5)

vln_save_to <- file.path(plot_dir, paste0("pbmc_cite.vlnplot_by_Deffner_etal_2024_makers.pdf"))
p <- (VlnPlot(pbmc_cite[, selected_cells], features = Deffner_etal_2024_makers, ncol = 4, pt.size = 0.1) &
  labs(x = NULL, y = NULL) &
  scale_x_discrete(limits = c("CD8 T", "CD4 T", "NK", "Mono", "B", "DC", "other T", "other")) &
  theme(title = element_text(size = 11), axis.text = element_text(size = 9))) + plot_layout(guides = "collect") 
ggsave(vln_save_to, plot = p, width = 8.25, height = 7.25)

dot_save_to <- file.path(plot_dir, paste0("pbmc_cite.dotplot_by_Deffner_etal_2024_makers.pdf"))
p <- DotPlot(pbmc_cite[, selected_cells], features = Deffner_etal_2024_makers) +
  scale_color_gradient(low = rb_gradient[5], high = rb_gradient[1]) +
  scale_x_discrete(limits = rev(Deffner_etal_2024_makers)) +
  scale_y_discrete(limits = c("CD8 T", "CD4 T", "NK", "Mono", "B", "DC", "other T", "other")) +
  labs(y = NULL) +
  RotatedAxis(0) +
  coord_flip()
ggsave(dot_save_to, plot = p, width = 5, height = 4)


# Find markers, either conserved or unique
DefaultAssay(pml_int) <- "SCT"
Idents(pml_int) <- "predicted.celltype.l1"

force <- TRUE
for (per_group in c("Responder", "Non-responder")) {
  selected_cells <- pml_int@meta.data %>%
    dplyr::filter(Timepoint == "BL", predicted.celltype.l1 == "CD8 T", Response == per_group | orig.ident == "Lymphoma_related_PML") %>%
    rownames()

  ceg_tbl_save_to <- file.path(object_dir, paste0("pml_integrated.", per_group, ".cd8_t.conserved_markers.txt"))
  if (!file.exists(ceg_tbl_save_to) || force) {
    ceg_tbl <- pml_int[, selected_cells] %>%
      FindConservedMarkers(ident.1 = "PMLCITEseq", ident.2 = "Lymphoma_related_PML", grouping.var = "predicted.celltype.l1", logfc.threshold = 0) %>%
      tibble::rownames_to_column("gene") %>%
      dplyr::mutate(group = paste0(per_group, ".PBMC_vs_CSF"))
    fwrite(ceg_tbl, ceg_tbl_save_to, sep = "\t")
  } else {
    ceg_tbl <- fread(ceg_tbl_save_to)
  }

  ego_save_to <- file.path(object_dir, paste0("pml_integrated.", per_group, ".cd8_t.ego.rds"))
  if (!file.exists(ego_save_to) || force) {
    selected_genes <- ceg_tbl %>%
      dplyr::filter(!stringr::str_detect(gene, "^RP[LS]|MT-"), !gene %in% chrom_xy_genes, `CD8 T_p_val_adj` < 0.05, abs(`CD8 T_avg_log2FC`) > 0.5) %>%
      dplyr::pull(gene)

    ego <- enrichGO(
      gene = selected_genes, universe = rownames(pml_int), keyType = "SYMBOL", OrgDb = "org.Hs.eg.db", ont = "ALL",
      pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05, readable = TRUE
    )
    saveRDS(ego, ego_save_to)
  } else {
    ego <- readRDS(ego_save_to)
  }
}

# Responders
ego_rs <- file.path(object_dir, "pml_integrated.Responder.cd8_t.ego.rds") %>% readRDS()
pdf(file.path(plot_dir, paste0("pml_integrated.Responder.cd8_t.ego.pdf")), width = 10, height = 10)
cnetplot(ego_rs, node_label = "category", showCategory = 15, layout = "fr", color.params = list(category='firebrick'))
dev.off()

# Non-responder
ego_nr <- file.path(object_dir, "pml_integrated.Non-responder.cd8_t.ego.rds") %>% readRDS()
pdf(file.path(plot_dir, paste0("pml_integrated.Non-responder.cd8_t.ego.pdf")), width = 10, height = 10)
cnetplot(ego_nr, node_label = "category", showCategory = 15, layout = "fr", color.params = list(category='firebrick'))
dev.off()

# Number of shared functional enrichment terms
go_term_venn <- dplyr::bind_rows(dplyr::mutate(ego_rs@result, group = "Responder"), dplyr::mutate(ego_nr@result, group = "Non-responder")) %>%
  dplyr::group_by(group) %>%
  dplyr::summarize(gene_list = list(ID)) %>%
  dplyr::pull(gene_list, group)
p <- ggVennDiagram(go_term_venn, label_size = 3) +
  scale_fill_gradient(low = rb_gradient[5], high = rb_gradient[1]) +
  coord_flip()
ggsave(file.path(plot_dir, "pml_integrated.conserved_go_term_venn.pdf"), plot = p, width = 7.25 * .6, height = 4.75 * .6)


shared_goid <- intersect(ego_nr@result$ID, ego_rs@result$ID)
# Responder specific
ego_rs@result$Description <- ego_rs@result$Description %>% purrr::map(~stringr::str_replace_all(.x, "regulation", "reg.")) %>% purrr::map(~stringr::str_replace_all(.x, "negative", "neg.")) %>% purrr::map(~stringr::str_replace_all(.x, "positive", "pos.")) %>% unlist()
selected_terms <- ego_rs@result %>% dplyr::filter(!ID %in% shared_goid, ONTOLOGY == "BP") %>% dplyr::arrange(p.adjust) %>% dplyr::pull(Description) %>% head(25)
pdf(file.path(plot_dir, paste0("pml_integrated.Responder_specific.cd8_t.ego.biological_processes.pdf")), width = 10, height = 10)
cnetplot(ego_rs, node_label = "category", showCategory = selected_terms, layout = "kk", color.params = list(category='firebrick'))
dev.off()

# Non-responder specific
ego_nr@result$Description <- ego_nr@result$Description %>% purrr::map(~stringr::str_replace_all(.x, "regulation", "reg.")) %>% purrr::map(~stringr::str_replace_all(.x, "negative", "neg.")) %>% purrr::map(~stringr::str_replace_all(.x, "positive", "pos.")) %>% unlist()
selected_terms <- ego_nr@result %>% dplyr::filter(!ID %in% shared_goid, ONTOLOGY == "BP") %>% dplyr::arrange(p.adjust) %>% dplyr::pull(Description) %>% head(25)
pdf(file.path(plot_dir, paste0("pml_integrated.Non-responder_specific.cd8_t.ego.biological_processes.pdf")), width = 10, height = 10)
cnetplot(ego_nr, node_label = "category", showCategory = selected_terms, layout = "kk", color.params = list(category='firebrick'))
dev.off()

# Shared
ego_nr@result$Description <- ego_nr@result$Description %>% purrr::map(~stringr::str_replace_all(.x, "regulation", "reg.")) %>% purrr::map(~stringr::str_replace_all(.x, "negative", "neg.")) %>% purrr::map(~stringr::str_replace_all(.x, "positive", "pos.")) %>% unlist()
selected_terms <- ego_nr@result %>% dplyr::filter(ID %in% shared_goid, ONTOLOGY == "BP") %>% dplyr::arrange(p.adjust) %>% dplyr::pull(Description) %>% head(25)
pdf(file.path(plot_dir, paste0("pml_integrated.shared_goid.cd8_t.ego.biological_processes.pdf")), width = 10, height = 10)
cnetplot(ego_nr, node_label = "category", showCategory = selected_terms, layout = "fr", color.params = list(category='firebrick'))
dev.off()



# Pseudo-time of base-line cells
selected_timepoint <- "BL"
selected_timepoint <- c("6W", "3M")
saving_flag <- ifelse(selected_timepoint %in% "BL", "base_line", "post_treatment") %>% unique()

cds_save_to <- file.path(object_dir, paste0("pml_integrated.cd8_t.", saving_flag, ".monocle3_cds"))
if (dir.exists(cds_save_to)) {
  pml_cds <- load_monocle_objects(cds_save_to)
} else {
  selected_cells <- pml_int@meta.data %>%
    dplyr::filter(Timepoint %in% selected_timepoint | orig.ident == "Lymphoma_related_PML", predicted.celltype.l1 == "CD8 T") %>% rownames()
  save_prefix <- file.path(object_dir, paste0("counts_cd8_t.", saving_flag), "outs", "filtered_feature_bc_matrix")
  write10xCounts(save_prefix, pml_int[, selected_cells]@assays$RNA@counts, version = "3", overwrite = TRUE)

  pml_cds <- load_cellranger_data(file.path(object_dir, paste0("counts_cd8_t.", saving_flag))) # Load 10X matrix

  # Add metadata
  colData(pml_cds) <- dplyr::left_join(
    as.data.frame(colData(pml_cds)),
    tibble::rownames_to_column(pml_int[, selected_cells]@meta.data, "CellBarcodes"),
    by = c("barcode" = "CellBarcodes")
  ) %>%
    tibble::column_to_rownames("barcode") %>%
    dplyr::mutate(barcode = rownames(.)) %>%
    dplyr::relocate(barcode) %>%
    as("DataFrame")

  pml_cds <- preprocess_cds(pml_cds) # Preprocess
  pml_cds <- align_cds(pml_cds, alignment_group = "SequencingPool") # Remove batch effects

  reducedDim(pml_cds, "PCA") <- pml_int[, selected_cells]@reductions$pca@cell.embeddings
  reducedDim(pml_cds, "UMAP") <- pml_int[, selected_cells]@reductions$umap@cell.embeddings

  pml_cds <- cluster_cells(pml_cds) # Cluster cells, based on UMAP
  pml_cds <- learn_graph(pml_cds) # Learn graph

  save_monocle_objects(pml_cds, file.path(object_dir, paste0("pml_integrated.cd8_t.", saving_flag, ".monocle3_cds")))
}

pml_cds$`Cell type` <- pml_cds$predicted.celltype.l2
pml_cds$`Dataset` <- lapply(pml_cds$orig.ident, function(x) {ifelse(x == "PMLCITEseq", "PBMC", "CSF")}) %>% unlist()

root_node <- ifelse(selected_timepoint %in% "BL", "Y_10", "Y_21") %>% unique()
pml_cds <- order_cells(pml_cds, root_pr_nodes = root_node)
p_timeseries <- plot_cells(pml_cds, color_cells_by = "pseudotime", label_leaves = FALSE, label_branch_points = FALSE) +
  scale_color_gradient(name = "Pseudo-time", low = "black", high = "yellow") +
  labs(subtitle = "Pseudotime", x = NULL, y = NULL) +
  lims(y = c(-2, NA), x = c(NA, 10)) +
  theme(legend.position = "right")

p_celltype <- plot_cells(pml_cds, color_cells_by = "Cell type", label_leaves = FALSE, label_branch_points = FALSE, label_cell_groups = FALSE) +
  scale_color_npg() +
  labs(subtitle = "Celltypes", x = NULL, y = NULL) +
  lims(y = c(-2, NA), x = c(NA, 10)) +
  theme(legend.position = "right")

dataset_label <- ifelse(selected_timepoint %in% "BL", "Dataset", "Timepoint")
p_dataset <- plot_cells(pml_cds, color_cells_by = dataset_label, label_leaves = FALSE, label_branch_points = FALSE, label_cell_groups = FALSE) +
  scale_color_jama() +
  labs(subtitle = dataset_label, x = NULL, y = NULL) +
  lims(y = c(-2, NA), x = c(NA, 10)) +
  theme(legend.position = "right")

selected_ps_features <- c("CCL4", "CCL5", "CCR7", "GNLY", "IFNG", "OASL", "LY9", "STXBP2")
p_gene_pdt <- plot_genes_in_pseudotime(pml_cds[selected_ps_features, ], color_cells_by = "predicted.celltype.l2", min_expr=0.1) +
  scale_color_npg() +
  labs(subtitle = "Marker gene", color = "Cell type") +
  theme(legend.position = "none")

p <- ((((p_dataset | p_celltype) / p_timeseries) + plot_layout(heights = c(1, 2))) | p_gene_pdt) + plot_layout(guides = "collect", widths = c(3, 1))
ggsave(file.path(plot_dir, paste0("pml_integrated.cd8_t.", saving_flag, ".pseudotime.pdf")), plot = p, width = 8, height = 6)
