#!/usr/bin/env Rscript
# File: atac_seq.r
# Author: Zhenhua Zhang
# E-mail: zhenhua.zhang217@gmail.com
# Created: Aug 17, 2023
# Updated: Aug 17, 2023

options(stringsAsFactors = FALSE, data.table.verbose = FALSE, error = rlang::last_trace)
suppressPackageStartupMessages({
  library(lobstr)
  library(data.table)
  library(tidyverse)
  library(patchwork)

  library(Seurat)
  library(Signac)
  library(Azimuth)
  library(SeuratData)
  library(SeuratDisk)
  library(CellMixS)
  library(SingleCellExperiment)
  library(harmony) # For integration

  library(GenomicRanges)
  library(EnsDb.Hsapiens.v86)
  library(BSgenome.Hsapiens.UCSC.hg38)
})


projdir <- "~/Documents/projects/wp_pml"
plot_dir <- file.path(projdir, "outputs/analysis/ATAC_seq/plots")
object_dir <- file.path(projdir, "outputs/analysis/ATAC_seq/objects")

overwrite <- FALSE
use_harmony <- TRUE
use_aggr_features <- TRUE

anchor_by <- "lsi" # rpca, cca, rlsi
normalize_by <- "TFIDF"

selected_pools <- c("3", "4")
pool_token <- paste(selected_pools, collapse = "_")
save_token <- paste0("pbmc.pool_", pool_token, ".archors_", anchor_by, ifelse(use_harmony, ".use_harmony", ""), ifelse(use_aggr_features, ".use_aggr_features", ""))

tar_cell_types <- c("Mono", "CD4 T", "CD8 T", "B", "NK", "DC")

# Shared markers
# TYROBP, (monocyte, NK cells), TRAC, (CD4+ T, CD8+ T), TMSB10, (CD4+ T, CD8+ T), CD3D, (CD4+ T, CD8+ T, Other T), CD3G, (CD4+ T, CD8+ T)
# NKG7, (NK, Other T), CST7, (NK, Other T), CD74, (DC, B), HLA-DQA1, (B, DC)
celltype_markers <- c(
  "CTSS", "FCN1", "LYZ", "PSAP", "S100A9", "AIF1", "SERPINA1", "CD14", "FCGR3A", # "TYROBP", "NEAT1", "MNDA", Monocytes
  "IL7R", "MAL", "LTB", "CD4", "LDHB", "CD3D", "CD3G", "TRAC", # "TMSB10", "TPT1", CD4+ T
  "CD8B", "CD8A", "HCST", "LINC02446", "CTSW", "CD3E", # "TMSB10", "CD3G", "CD3D", "TRAC", CD8+ T
  "TRDC", "GZMK", "KLRB1", "TRGC2", "LYAR", "KLRG1", "GZMA", # "CD3D", "CST7", "NKG7", Other T
  "KLRD1", "GNLY", "PRF1", "CD247", "KLRF1", "GZMB", # "TYROBP", "FCER1G", "CST7", NK
  "CD79A", "CD79B", "RALGPS2", "MS4A1", "BANK1", "TNFRSF13C", "IGHM", "MEF2C", # "CD74", "HLA-DQA1", B
  "HLA-DPA1", "HLA-DPB1", "CCDC88A", "HLA-DRA", "HLA-DMA", "CST3", "HLA-DQB1", "HLA-DRB1" # "CD74", "HLA-DQA1", DC
)


# Annotations
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- "UCSC"
genome(annotations) <- "hg38"

# Load ATAC-seq data
pools <- paste0("ATACB2P", selected_pools) %>% rlang::set_names(.)
bcmat_path <- file.path(projdir, "outputs/readcounts/ATAC_seq", pools, "outs") %>% rlang::set_names(pools)
dmres_path <- file.path(projdir, "outputs/demultiplex/ATAC_seq", pools, "vireo_outdir_ref/donor_ids.tsv") %>% purrr::set_names(pools)


# Loading data using features from the pool per se
if (FALSE) { # Quick test
  pbmc_list_raw <- lapply(bcmat_path, function(pp) {
    pn <- stringr::str_split(pp, pattern = "/", simplify = TRUE) %>% purrr::keep(~stringr::str_detect(.x, "ATACB2P"))
    cat("[I]: Creating from raw data ...\n")
    metadata <- read.csv(file = file.path(pp, "singlecell.csv"), header = TRUE, row.names = 1) %>% dplyr::mutate(pool = pn)
    counts <- Read10X_h5(filename = file.path(pp, "filtered_peak_bc_matrix.h5"))
    po <- CreateChromatinAssay(counts = counts, sep = c(":", "-"), fragments = file.path(pp, "fragments.tsv.gz"), min.cells = 10, min.features = 200) %>% CreateSeuratObject(assay = "ATAC", meta.data = metadata)

    Annotation(po) <- annotations # Add annotations
    po <- TSSEnrichment(object = po, fast = FALSE) # compute TSS enrichment score per cell
    po <- NucleosomeSignal(object = po) # compute nucleosome signal score per cell

    return(invisible(po))
  })

  pbmc_raw <- lapply(pools, function(pp) { pbmc_list_raw[[pp]] %>% RunTFIDF() %>% FindTopFeatures(min.cutoff = 10) %>% RunSVD() })
  pbmc_raw_merge <- merge(pbmc_raw[[1]], pbmc_raw[2:length(pbmc_raw)]) %>% FindTopFeatures(min.cutoff = 500) %>% RunTFIDF() %>% RunSVD() %>% RunUMAP(reduction = "lsi", dims = 2:20) 
  pbmc_raw_int <- RunHarmony(pbmc_raw_merge, c("pool"), dims.use = 2:30, reduction = "lsi", assay.use = "ATAC", project.dim = FALSE) # Run HARMONY
  pbmc_raw_int <- RunUMAP(pbmc_raw_int, reduction = "harmony", dims = 2:30)
  save_to <- file.path(plot_dir, "integrated", paste0("pbmc.raw.integrated.umap.by_pool.pdf"))
  p <- DimPlot(pbmc_raw_int, reduction = "umap", group.by = "pool")
  ggsave(save_to, plot = p)
}


# Create common features
all_chroms <- paste0("chr", c(1:22, "X", "Y"))
if (use_aggr_features) {
# Create common features from `cellranger-atac aggr` pipeline
  common_features <- file.path(projdir, "outputs/aggregate", paste0("ATACseq_pool_", pool_token), "outs/peaks.bed") %>%
    read.table(col.names = c("chr", "start", "end")) %>%
    dplyr::filter(chr %in% all_chroms) %>%
    makeGRangesFromDataFrame() %>%
    (function(gr) { w <- width(gr); gr[w < 2000 & w > 100]})
} else {
# Create common features from peaking calling results of each pool
  common_features <- lapply(bcmat_path, function(pp, .overwrite = FALSE) {
    pn <- stringr::str_split(pp, pattern = "/", simplify = TRUE) %>% purrr::keep(~stringr::str_detect(.x, "ATACB2P") & .x %in% pools)
    peak_gr <- file.path(pp, "peaks.bed") %>%
      read.table(col.names = c("chr", "start", "end")) %>%
      dplyr::filter(chr %in% all_chroms) %>%
      makeGRangesFromDataFrame()
  }) %>%
    Reduce(c, .) %>%
    IRanges::reduce(ignore.strand = TRUE) %>%
    (function(gr) { w <- width(gr); gr[w < 2000 & w > 100]})
}


# Create fratments objects, using common features
pbmc_list <- lapply(pools, function(pn, .dmres_path, .bcmat_path, .common_features, .overwrite) {
  save_to <- file.path(object_dir, "raw", paste0("pbmc.", pn, ".raw.rds"))
  if (!file.exists(save_to) || .overwrite)  {
    cat("[I]: Creating from raw data ...\n")
    dm_tab <- fread(.dmres_path[pn]) %>%
      dplyr::select(
        Cellbarcode = cell, Vireo_assignment = donor_id, Vireo_prob_max = prob_max, Vireo_prob_doublet = prob_doublet,
        Vireo_dbl_llr = doublet_logLikRatio
      )
    metadata <- read.csv(file.path(.bcmat_path[pn], "filtered_peak_bc_matrix", "barcodes.tsv"), col.names = "Cellbarcode") %>%
      dplyr::mutate(pool = pn) %>%
      dplyr::left_join(dm_tab, by = "Cellbarcode") %>%
      tibble::column_to_rownames("Cellbarcode")

    frag_counts <- CountFragments(fragments = file.path(.bcmat_path[pn], "fragments.tsv.gz"))
    atac_cells <- frag_counts[frag_counts$frequency_count >= 1000, "CB"]
    fragments <- CreateFragmentObject(path = file.path(.bcmat_path[pn], "fragments.tsv.gz"), cells = atac_cells)
    counts <- FeatureMatrix(fragments = fragments, features = .common_features, cells = atac_cells)
    po <- CreateChromatinAssay(counts = counts, fragments = fragments, min.features = 500, min.cells = 5) %>%
      CreateSeuratObject(assay = "ATAC", meta.data = metadata)

    Annotation(po) <- annotation # Add annotations
    po <- TSSEnrichment(object = po, fast = FALSE) # compute TSS enrichment score per cell
    po <- NucleosomeSignal(object = po) # compute nucleosome signal score per cell

    p <- DensityScatter(po, x = 'nCount_ATAC', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE) # Density scatter plot to show
    file.path(plot_dir, "quality_control", paste0(pn, ".pre_qc.density_scatter.tss_enrichment.pdf")) %>% ggsave(plot = p)

    p <- DensityScatter(po, x = 'nCount_ATAC', y = 'nucleosome_signal', log_x = TRUE, quantiles = TRUE) # Density scatter plot to show
    file.path(plot_dir, "quality_control", paste0(pn, ".pre_qc.density_scatter.nucleosome_signal.pdf")) %>% ggsave(plot = p)

    p <- DensityScatter(po, x = 'nCount_ATAC', y = 'nFeature_ATAC', log_x = TRUE, quantiles = TRUE) # Density scatter plot to show
    file.path(plot_dir, "quality_control", paste0(pn, ".pre_qc.density_scatter.n_feature_atac.pdf")) %>% ggsave(plot = p)

    po$high.tss <- ifelse(po$TSS.enrichment > 1, 'High', 'Low') # TSS enrichment, using 3 as threshold
    p <- TSSPlot(po, group.by = 'high.tss') + NoLegend()
    file.path(plot_dir, "quality_control", paste0(pn, ".pre_qc.tss_plot.pdf")) %>% ggsave(plot = p)

    po$nucleosome_group <- ifelse(po$nucleosome_signal > 4, 'NS > 4', 'NS < 4') # cells grouped by high or low nucleosomal signal strength, using 2 as threshold
    p <- FragmentHistogram(object = po, group.by = 'nucleosome_group')
    file.path(plot_dir, "quality_control", paste0(pn, ".pre_qc.fragment_histogram.pdf")) %>% ggsave(plot = p)

    p <- VlnPlot(object = po, features = c('nCount_ATAC', "nFeature_ATAC", 'TSS.enrichment', 'nucleosome_signal', pt.size = 0.1, ncol = 5)) # QC metric separately using a violin plot
    file.path(plot_dir, "quality_control", paste0(pn, ".pre_qc.violin_plot.pdf")) %>% ggsave(plot = p)

    cat("[I]: Dumping into disk ...\n"); saveRDS(po, save_to)
  } else {
    cat("[I]: Loading from disk ...\n"); po <- readRDS(save_to)
  }

  return(invisible(po))
}, .bcmat_path = bcmat_path, .dmres_path = dmres_path, .common_features = common_features, .overwrite = FALSE)


# Quality control figures in one plot.
save_to <- file.path(plot_dir, "quality_control", paste0("pbmc.pool_", pool_token, ".control_metrics.violin.pdf"))
if (!file.exists(save_to) || overwrite) {
  meta_data_all <- lapply(pools, function(pp) { as.data.frame(pbmc_list[[pp]]@meta.data) %>% dplyr::mutate(pool = pp, barcode = rownames(.)) }) %>%
    Reduce(rbind, .) %>%
    dplyr::mutate(cellbarcode = rownames(.)) %>%
    dplyr::select(cellbarcode, pool, nCount_ATAC, nFeature_ATAC, TSS.enrichment, nucleosome_signal) %>%
    tidyr::pivot_longer(cols = c(nCount_ATAC, nFeature_ATAC, TSS.enrichment, nucleosome_signal)) %>%
    dplyr::mutate(name = factor(name, levels = c("nCount_ATAC", "nFeature_ATAC", "TSS.enrichment", "nucleosome_signal")))

  p <- ggplot(data = meta_data_all) +
    geom_point(aes(x = pool, y = value), position = "jitter", size = 0.05) +
    geom_violin(aes(x = pool, y = value, fill = pool), color = "white", alpha = 0.75, linewidth = 0.2) +
    facet_wrap(~name, scales = "free_y", nrow = 1) +
    scale_fill_discrete(name = "", type = c("darkred", "red", "tan3", "orange")) +
    theme_classic() +
    theme(axis.text.x.bottom = element_text(angle = 30, hjust = 1, vjust = 1), strip.background = element_rect(linetype = "blank"), legend.position = "top") +
    labs(x = NULL, y = NULL)

  if (!file.exists(save_to) || overwrite) {
    ggsave(save_to, plot = p, width = 9, height = 4)
  }
}


# QC and normalization
# TODO: add other metadata, e.g. age, gender, and disease group, etc.
hard_threshold <- tibble::tribble( # Hard threshold per pool
  ~pool, ~min_ncount_atac, ~max_ncount_atac, ~min_nfeature_atac, ~max_nfeature_atac, ~min_tss_enrichment, ~max_tss_enrichment, ~max_nucleosome_signal,
  # "ATACB2P1", 500, 2000, 500, 5000, 0.75, 2, 4,
  # "ATACB2P2", 500, 2000, 750, 2000, 2.5, 9, 4,
  "ATACB2P3", 1223, 9349, 1010, 9045, 2.83, 6.67, 2.09,
  "ATACB2P4", 1143, 5408, 700, 4598, 1.87, 5.71, 2.98,
)

save_to <- file.path(object_dir, "post_qc", paste0("pbmc.pool_", pool_token, ".norm_", normalize_by, ".rds"))
if (!file.exists(save_to) || overwrite) {
  pbmc <- lapply(pools, function(pn, .pbmc, .threshold) {
    po <- .pbmc[[pn]]
    po$pool <- pn
    po$dataset <- pn

    params <- .threshold %>% dplyr::filter(pool == pn)
    min_ncount_atac <- params["min_ncount_atac"] %>% as.integer()
    max_ncount_atac <- params["max_ncount_atac"] %>% as.integer()
    min_nfeature_atac <- params["min_nfeature_atac"] %>% as.integer()
    max_nfeature_atac <- params["max_nfeature_atac"] %>% as.integer()
    min_tss_enrichment <- params["min_tss_enrichment"] %>% as.double()
    max_tss_enrichment <- params["max_tss_enrichment"] %>% as.double()
    max_nucleosome_signal <- params["max_nucleosome_signal"] %>% as.double()
    tar_cells <- as.data.frame(po@meta.data) %>%
      dplyr::filter(min_ncount_atac <= nCount_ATAC, nCount_ATAC <= max_ncount_atac, min_nfeature_atac <= nFeature_ATAC, nFeature_ATAC <= max_nfeature_atac, min_tss_enrichment <= TSS.enrichment,  TSS.enrichment <= max_tss_enrichment, nucleosome_signal <= max_nucleosome_signal) %>%
      dplyr::filter(!Vireo_assignment %in% c("doublet", "unassigned"), !is.na(Vireo_assignment)) %>%
      rownames()

    p <- DensityScatter(po, x = 'nCount_ATAC', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE) # Density scatter plot to show
    file.path(plot_dir, "quality_control", paste0(pn, ".post_qc.density_scatter.tss_enrichment.pdf")) %>% ggsave(plot = p)

    p <- DensityScatter(po, x = 'nCount_ATAC', y = 'nucleosome_signal', log_x = TRUE, quantiles = TRUE) # Density scatter plot to show
    file.path(plot_dir, "quality_control", paste0(pn, ".post_qc.density_scatter.nucleosome_signal.pdf")) %>% ggsave(plot = p)

    p <- DensityScatter(po, x = "nCount_ATAC", y = 'nFeature_ATAC', log_x = TRUE, quantiles = TRUE) # Density scatter plot to show
    file.path(plot_dir, "quality_control", paste0(pn, ".post_qc.density_scatter.n_feature_atac.pdf")) %>% ggsave(plot = p)

    p <- TSSPlot(po, group.by = "high.tss") + NoLegend()
    file.path(plot_dir, "quality_control", paste0(pn, ".post_qc.tss_plot.pdf")) %>% ggsave(plot = p)

    p <- FragmentHistogram(object = po, group.by = "nucleosome_group")
    file.path(plot_dir, "quality_control", paste0(pn, ".post_qc.fragment_histogram.pdf")) %>% ggsave(plot = p)

    p <- VlnPlot(object = po, features = c('nCount_ATAC', 'TSS.enrichment', 'nucleosome_signal', pt.size = 0.1, ncol = 5)) # QC metric separately using a violin plot
    file.path(plot_dir, "quality_control", paste0(pn, ".post_qc.violin_plot.pdf")) %>% ggsave(plot = p)

    po <- po[, tar_cells] %>% RunTFIDF(verbose = FALSE) %>% FindTopFeatures(min.cutoff = 10) %>% RunSVD(verbose = FALSE)

    p <- DepthCor(po) # Correlation between each LSI component and sequencing depth. To remove LSI components highly corrlated with sequencing depth
    file.path(plot_dir, "quality_control", paste0(pn, ".depth_cor.pdf")) %>% ggsave(plot = p)

    return(invisible(po))
  }, .pbmc = pbmc_list, .threshold = hard_threshold)

  cat("[I]: Dumping into disk ...\n"); saveRDS(pbmc, save_to)
} else {
  cat("[I]: Loading from disk ...\n"); pbmc <- readRDS(save_to)
}


#
## Integration
#
if (length(pools) > 1) {
# Merge to obtain a reduction
  save_to <- file.path(object_dir, "merged", paste0("pbmc.pool_", pool_token, ".merged.rds"))
  if (!file.exists(save_to) || overwrite) {
    cat("[I]: Merging objects ...\n")
    pbmc_merge <- merge(pbmc[[1]], pbmc[2:length(pbmc)]) %>% RunTFIDF(verbose = FALSE) %>% FindTopFeatures(min.cutoff = 50) %>% RunSVD() %>% RunUMAP(reduction = "lsi", dims = 2:20) 
    cat("[I]: Dumping into disk ...\n"); saveRDS(pbmc_merge, save_to)
  } else {
    cat("[I]: Loading from disk ...\n"); pbmc_merge <- readRDS(save_to)
  }

  # Integration
  save_to <- file.path(object_dir, "integrated", paste0(save_token, ".integrated.rds"))
  if (!file.exists(save_to) || overwrite) {
    if (use_harmony) {
      cat("[I]: Integrating objects using HARMONY ...\n")
      pbmc_int <- RunHarmony(pbmc_merge, c("pool"), dims.use = 2:30, reduction = "lsi", assay.use = "ATAC", project.dim = FALSE) # Run HARMONY
    } else {
      cat("[I]: Integrating objects using Seurat internal method ...\n")
      int_anchors <- FindIntegrationAnchors(pbmc, anchor.features = VariableFeatures(pbmc_merge), reduction = "rlsi", dims = 2:30)
      pbmc_int <- IntegrateEmbeddings(int_anchors, reductions = pbmc_merge[["lsi"]], k.weight = 50, new.reduction.name = "lsi", dims.to.integrate = 1:30)
    }
    cat("[I]: Dumping into disk ...\n"); saveRDS(pbmc_int, save_to)
  } else {
    cat("[I]: Loading from disk ...\n"); pbmc_int <- readRDS(save_to)
  }
} else {
  pbmc_int <- pbmc[[1]]
}


# Reduction by UMAP
umap_reduction <- ifelse(use_harmony, "harmony", "lsi")
pbmc_int <- RunUMAP(pbmc_int, reduction = umap_reduction, dims = 2:10)
save_to <- file.path(plot_dir, "integrated", paste0(save_token, ".umap_by_pool.pdf"))
if (!file.exists(save_to) || overwrite) {
  p <- DimPlot(pbmc_int, reduction = "umap", group.by = "pool"); ggsave(save_to, plot = p)
}

# Identify clusters
pbmc_int <- FindNeighbors(pbmc_int, reduction = umap_reduction, dims = 2:10)
pbmc_int <- FindClusters(pbmc_int, algorithm = 3)

save_to <- file.path(plot_dir, "integrated", paste0(save_token, ".by_seruat_clusters.pdf"))
if (!file.exists(save_to) || overwrite) {
  p <- DimPlot(pbmc_int, reduction = "umap", group.by = "seurat_clusters")
  ggsave(save_to, plot = p)
}

save_to <- file.path(object_dir, "integrated", paste0(save_token, ".clustered.rds"))
if (!file.exists(save_to) || overwrite) {
  cat("[I]: Dumping into disk ...\n"); saveRDS(pbmc_int, save_to)
} else {
  cat("[I]: Loading from disk ...\n"); pbmc_int <- readRDS(save_to)
}


# Add the gene activity matrix to the Seurat object as a new assay and normalize it
gene_activities <- GeneActivity(pbmc_int)
pbmc_int[['RNA']] <- CreateAssayObject(counts = gene_activities)
pbmc_int <- NormalizeData(pbmc_int, assay = 'RNA', normalization.method = 'LogNormalize', scale.factor = median(pbmc_int$nCount_RNA))

save_to <- file.path(plot_dir, "integrated", paste0(save_token, ".integrated.celltype_markers.dot_plot.pdf"))
if (!file.exists(save_to) || overwrite) {
  DefaultAssay(pbmc_int) <- "RNA"
  p <- DotPlot(pbmc_int, features = celltype_markers) + scale_size(range = c(0.2, 6.5)) + coord_flip() + RotatedAxis()
  ggsave(save_to, p, width = 6, height = 10)
}

save_to <- file.path(object_dir, "integrated", paste0(save_token, ".clustered.annotated.rds"))
if (!file.exists(save_to) || overwrite) {
  cat("[I]: Dumping into disk ...\n"); saveRDS(pbmc_int, save_to)
} else {
  cat("[I]: Loading from disk ...\n"); pbmc_int <- readRDS(save_to)
}





# -------------------------- To be removed -------------------------------
# Per pool
for (pid in as.character(1:4)) {
  save_to <- file.path(object_dir, "per_pool", paste0("pbmc.pool_", tar_pool, ".rds"))
  if (!file.exists(save_to) || overwrite) {
    tar_pool <- paste0("ATACB2P", pid) %>% rlang::set_names(.)
    bcmat_path <- file.path(projdir, "outputs/readcounts/ATAC_seq", tar_pool, "outs") %>% rlang::set_names(tar_pool)

    counts <- Read10X_h5(filename = file.path(bcmat_path, "filtered_peak_bc_matrix.h5"))
    metadata <- read.csv(file = file.path(bcmat_path, "singlecell.csv"), header = TRUE, row.names = 1)
    chrom_assay <- CreateChromatinAssay(counts = counts, sep = c(":", "-"), fragments = file.path(bcmat_path, "fragments.tsv.gz"), min.cells = 10, min.features = 200)
    pbmc_pp <- CreateSeuratObject(counts = chrom_assay, assay = "peaks", meta.data = metadata)
    Annotation(pbmc_pp) <- annotation

    pbmc_pp <- NucleosomeSignal(object = pbmc_pp)
    pbmc_pp <- TSSEnrichment(object = pbmc_pp, fast = FALSE)

    params <- hard_threshold %>% dplyr::filter(pool == tar_pool)
    min_ncount_atac <- params["min_ncount_atac"] %>% as.integer()
    max_ncount_atac <- params["max_ncount_atac"] %>% as.integer()
    min_tss_enrichment <- params["min_tss_enrichment"] %>% as.double()
    max_tss_enrichment <- params["max_tss_enrichment"] %>% as.double()
    max_nucleosome_signal <- params["max_nucleosome_signal"] %>% as.double()
    tar_cells <- pbmc_pp@meta.data %>%
      as.data.frame() %>%
      dplyr::filter(min_ncount_atac <= nCount_peaks, nCount_peaks <= max_ncount_atac, min_tss_enrichment <= TSS.enrichment,  TSS.enrichment <= max_tss_enrichment, nucleosome_signal <= max_nucleosome_signal) %>%
      rownames()

    pbmc_pp <- FindTopFeatures(pbmc_pp[, tar_cells], min.cutoff = 'q0')
    pbmc_pp <- RunTFIDF(pbmc_pp)
    pbmc_pp <- RunSVD(pbmc_pp)
    pbmc_pp <- RunUMAP(pbmc_pp, reduction = 'lsi', dims = 2:30)
    pbmc_pp <- FindNeighbors(pbmc_pp, reduction = 'lsi', dims = 2:30)
    pbmc_pp <- FindClusters(pbmc_pp, verbose = FALSE, algorithm = 3)

    fig_save_to <- file.path(plot_dir, "per_pool", paste0("pbmc.pool_", tar_pool, ".umap.pdf"))
    p <- DimPlot(pbmc_pp, label = TRUE) + NoLegend()
    ggsave(fig_save_to, plot = p, width = 6, height = 6)

    saveRDS(pbmc_pp, save_to)
  } else {
    cat("[I]: Found from disk ...\n")
  }
}
