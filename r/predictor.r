#!/usr/bin/env Rscript
# File: predictor.r
# Author: Zhenhua Zhang
# E-mail: zhenhua.zhang217@gmail.com
# Created: Nov 09, 2025
# Updated:

# Machine learning model to predict PML outcomes
options(stringsAsFactors = FALSE, data.table.verbose = FALSE, error = rlang::last_trace, future.globals.maxSize = 1024 * 10 ^ 10)
suppressPackageStartupMessages({
  library(lobstr)
  library(data.table)
  library(tidyverse)
  library(patchwork)

  library(kernelshap)
  library(caret)
  library(doParallel)
  library(MLeval)
  library(gbm)
  library(harmony)

  library(Seurat)
  library(SeuratData)
  library(SeuratDisk)

  library(ComplexHeatmap)
  library(clusterProfiler)
  library(circlize)
  library(RColorBrewer)
  library(org.Hs.eg.db)
})


meta_data <- function(o, as.df = TRUE) { if (as.df) as.data.frame(o@meta.data) else o@meta.data }

red_blue <- colorRampPalette(brewer.pal(11, "RdBu"))(11)
blues <- colorRampPalette(brewer.pal(9, "Blues"))(9)
col_fun <- colorRamp2(c(0, 20), c("#90D5FFBB", "darkblue"))


#
## Define parameters
#
projdir <- "~/Documents/projects/wp_pml"
plot_dir <- file.path(projdir, "outputs/analysis/prediction/plots")
object_dir <- file.path(projdir, "outputs/analysis/prediction/models")

# Remove genes on X and Y chromosomes.
chrom_xy_genes <- fread("~/Documents/projects/wp_pml/inputs/references/sextual_chromosome_genes.txt", header = F) %>% pull(V2)

all_donors <- c("PML0002", "PML0009", "PML0017", "PML0022", "PML0025")

#
## Load datasets
#
pbmc_cite <- "outputs/analysis/integrated/objects/pml_citeseq.rds" %>% file.path(projdir, .) %>% readRDS()
pbmc_rna <- "outputs/analysis/integrated/objects/pml_rnaseq.rds" %>% file.path(projdir, .) %>% readRDS() %>%
  subset(., subset = nCount_RNA >= 1000 & nCount_RNA <= 17000 & nFeature_RNA >= 600 & nFeature_RNA <= 4500)

# Removing batch effects
var_to_regress <- c("percent_mt", "S.Score", "G2M.Score")
pbmc_list <- c(SplitObject(pbmc_cite, "SequencingPool"), SplitObject(pbmc_rna, "SequencingPool")) %>%
  lapply(function(o) {
    DefaultAssay(o) <- "RNA"
    o[["SCT"]] <- NULL
    SCTransform(o, assay = "RNA", vst.flavor = "v2", vars.to.regress = var_to_regress, method = "glmGamPoi") %>%
      RunPCA(assay = "SCT", verbose = FALSE)
  })

integration_features <- SelectIntegrationFeatures(pbmc_list)
pbmc <- PrepSCTIntegration(pbmc_list, anchor.features = integration_features, verbose = FALSE)
pbmc_int <- FindIntegrationAnchors(
  pbmc, normalization.method = "SCT", anchor.features = integration_features, reduction = "rpca", dims = 1:30,
  k.anchor = 20, verbose = FALSE
)
pbmc_int <- IntegrateData(pbmc_int, normalization.method = "SCT", dims = 1:30, verbose = FALSE)
pbmc_int <- RunPCA(pbmc_int, assay = "integrated", verbose = FALSE)
pbmc_int <- RunUMAP(pbmc_int, assay = "integrated", dims = 1:30, verbose = FALSE)
pbmc_int <- FindNeighbors(pbmc_int, assay = "integrated", dims = 1:30, verbose = FALSE)
pbmc_int <- FindClusters(pbmc_int, graph.name = "integrated_nn", resolution = 0.1, verbose = FALSE)
Idents(pbmc_int) <- "predicted.celltype.l1"
pbmc_int <- PrepSCTFindMarkers(pbmc_int, verbose = FALSE)
rm(list = c("pbmc_list", "pbmc_cite", "pbmc_rna")); gc()


#
## Training
#
train_val_loop <- function(obj, train_cells, val_cells, saving_flag, object_dir) {
  pbmc_train <- subset(obj, cells = cells_train)
  pbmc_val <- subset(obj, cells = cells_val)

  #
  ## Train and evaluate the model
  #
  model_path <- file.path(object_dir, paste0("train.", saving_flag, ".gbm_fit.model.rds"))
  feature_list <- file.path(object_dir, paste0("train.", saving_flag, ".gbm_fit.feature_list.txt"))
  varimp_path <- file.path(object_dir, paste0("train.", saving_flag, ".gbm_eval.variable_importance.rds"))
  train_eval_metric_path <- file.path(object_dir, paste0("train.", saving_flag, ".gbm_eval.eval_metric.rds"))
  confusion_matrix_path <- file.path(object_dir, paste0("train.", saving_flag, ".gbm_eval.confusion_matrix.rds"))
  if (all(file.exists(model_path, feature_list, varimp_path, confusion_matrix_path))) {
    cat("[I]: Loading from disk ...\n")
    gbm_fit <- readRDS(model_path)
    selected_features <- fread(feature_list) %>% pull(features)
    test_confusion_matrix <- readRDS(confusion_matrix_path)
  } else {
    cat("[I]: Training model ...\n")
    # Define features and cells to work on
    selected_cells <- meta_data(pbmc_train) %>% dplyr::filter(Timepoint == "BL")
    if (which_cell_type == "all") {
      selected_cells <- selected_cells %>% rownames()
    } else {
      selected_cells <- selected_cells %>% dplyr::filter(predicted.celltype.l1 == which_cell_type) %>% rownames()
    }

    selected_features <- (rowSums(pbmc_train[, selected_cells]@assays$SCT@counts >= 3) >= 200) %>%
      purrr::keep(~.x) %>% names() %>% purrr::discard(~.x %in% chrom_xy_genes)
    # selected_features <- VariableFeatures(pbmc_train) %>% purrr::discard(~.x %in% chrom_xy_genes)
    data.frame(features = selected_features) %>% write.table(file = feature_list, row.names = F, quote = F)

    # Obtain whole dataset
    train_label_tbl <- pbmc_train[selected_features, selected_cells] %>% meta_data() %>%
      rownames_to_column("cell_barcode") %>% dplyr::select(cell_barcode, Response)
    train_data_tbl <- pbmc_train[selected_features, selected_cells]@assays$SCT@data %>%
      as.matrix() %>% t() %>% as.data.frame() %>% rownames_to_column("cell_barcode") %>%
      dplyr::left_join(train_label_tbl, by = "cell_barcode") %>%
      dplyr::mutate(Response = dplyr::if_else(Response == "Responder", "Responder", "Non_responder")) %>%
      dplyr::select(-cell_barcode)

    # Split dataset into training and testing
    in_training <- createDataPartition(train_data_tbl$Response, p = .75, list = FALSE)
    internal_train_tbl <- train_data_tbl[ in_training,]

    # Control the training parameters
    fit_control <- trainControl(method = "repeatedcv", number = 10, repeats = 10, classProbs = TRUE, summaryFunction = twoClassSummary)

    # Use multiple cores to train the model
    cl <- makePSOCKcluster(8)
    registerDoParallel(cl)
    gbm_fit <- train(Response ~ ., internal_train_tbl, method = "gbm", trControl = fit_control, metric = "ROC")
    stopCluster(cl)

    # Define the variable importance
    test_var_imp <- varImp(gbm_fit, scale = FALSE)
    saveRDS(test_var_imp, varimp_path)

    # Split dataset into training and testing
    cat("[I]: Calibrating model using important features ...\n")
    selected_cols <- as.data.frame(test_var_imp$importance) %>%
      tibble::rownames_to_column("Feature") %>%
      dplyr::filter(Overall <= 0) %>%
      dplyr::pull(Feature) %>%
      (function(x) which(! colnames(train_data_tbl) %in% x))
    internal_train_tbl <- train_data_tbl[ in_training, selected_cols]

    cl <- makePSOCKcluster(8)
    registerDoParallel(cl)
    gbm_fit <- train(Response ~ ., internal_train_tbl, method = "gbm", trControl = fit_control, metric = "ROC")
    stopCluster(cl)

    # Save the model
    saveRDS(gbm_fit, model_path) # Save the model

    # Evaluate model performance
    cat("[I]: Evaluating model performance ...\n")
    internal_test_tbl  <- train_data_tbl[-in_training, selected_cols]
    test_pred_tbl <- predict(gbm_fit, newdata = internal_test_tbl, type = "prob") %>%
      cbind(internal_test_tbl) %>%
      dplyr::mutate(pred = dplyr::case_when( Responder > 0.50 ~ "Responder", Responder <= 0.50 ~ "Non_responder")) %>%
      dplyr::select(obs = Response, pred, Y = Responder, N = Non_responder) %>%
      dplyr::mutate(obs = factor(obs, levels = c("Responder", "Non_responder"))) %>%
      dplyr::mutate(pred = factor(pred, levels = c("Responder", "Non_responder")))

    # Plot ROC
    test_eval_tbl <- dplyr::select(test_pred_tbl, Responder = Y, Non_responder = N, obs) %>%
      evalm(positive = "Responder")
    file.path(plot_dir, paste0("train.", saving_flag, ".gbm_eval.roc_curve.pdf")) %>%
      ggsave(., plot = test_eval_tbl[[1]], width = 8, height = 5)
    file.path(plot_dir, paste0("train.", saving_flag, ".gbm_eval.proc_curve.pdf")) %>%
      ggsave(., plot = test_eval_tbl[[2]], width = 8, height = 5)
    file.path(plot_dir, paste0("train.", saving_flag, ".gbm_eval.prg_curve.pdf")) %>%
      ggsave(., plot = test_eval_tbl[[3]], width = 8, height = 5)

    saveRDS(test_eval_tbl, train_eval_metric_path)

    # Confusion matrix
    test_confusion_matrix <- confusionMatrix(
      data = test_pred_tbl$pred, reference = test_pred_tbl$obs, mode = "everything", positive = "Responder"
    )
    saveRDS(test_confusion_matrix, confusion_matrix_path)
  }

  #
  ## Validate the model using external dataset.
  #
  # Preprocessing
  cat("[I]: Validating model ...\n")
  val_selected_cells <- meta_data(pbmc_val) %>% dplyr::filter(Timepoint == "BL")
  if (which_cell_type == "all") {
    val_selected_cells <- val_selected_cells %>% rownames()
  } else {
    val_selected_cells <- val_selected_cells %>% dplyr::filter(predicted.celltype.l1 == which_cell_type) %>% rownames()
  }

  val_response_tbl <- pbmc_val[, val_selected_cells] %>% meta_data() %>%
    rownames_to_column("cell_barcode") %>%
    dplyr::select(cell_barcode, Response) %>%
    dplyr::mutate(Response = dplyr::if_else(Response == "Responder", "Responder", "Non_responder"))

  val_data_tbl <- pbmc_val[selected_features, val_selected_cells]@assays$SCT@data %>%
    as.matrix() %>% t() %>% as.data.frame() %>% rownames_to_column("cell_barcode") %>%
    dplyr::left_join(val_response_tbl, by = "cell_barcode") %>%
    dplyr::select(-cell_barcode)

  val_pred_tbl <- predict(gbm_fit, newdata = val_data_tbl, type = "prob") %>%
    cbind(val_data_tbl) %>%
    dplyr::mutate(pred = dplyr::case_when(Responder > 0.50 ~ "Responder", Responder <= 0.50 ~ "Non_responder")) %>%
    dplyr::select(obs = Response, pred, Y = Responder, N = Non_responder) %>%
    dplyr::mutate(obs = factor(obs, levels = c("Responder", "Non_responder"))) %>%
    dplyr::mutate(pred = factor(pred, levels = c("Responder", "Non_responder")))

  # Validation performance
  val_confusion_matrix_path <- file.path(object_dir, paste0("validation.", saving_flag, ".gbm_eval.confusion_matrix.rds"))
  val_eval_metric_path <- file.path(object_dir, paste0("validation.", saving_flag, ".gbm_eval.eval_metric.rds"))
  if (all(file.exists(c(val_confusion_matrix_path, val_eval_metric_path)))) {
    val_confusion_matrix <- readRDS(val_confusion_matrix_path)
    val_eval_metric <- readRDS(val_eval_metric_path)
  } else {
    # Plotting validation results
    val_eval_tbl <- dplyr::select(val_pred_tbl, Responder = Y, Non_responder = N, obs) %>% evalm(positive = "Responder")
    file.path(plot_dir, paste0("validation.", saving_flag, ".gbm_eval.roc_curve.pdf")) %>%
      ggsave(., plot = val_eval_tbl[[1]], width = 8, height = 5)
    file.path(plot_dir, paste0("validation.", saving_flag, ".gbm_eval.proc_curve.pdf")) %>%
      ggsave(., plot = val_eval_tbl[[2]], width = 8, height = 5)
    file.path(plot_dir, paste0("validation.", saving_flag, ".gbm_eval.prg_curve.pdf")) %>%
      ggsave(., plot = val_eval_tbl[[3]], width = 8, height = 5)

    saveRDS(val_eval_tbl, val_eval_metric_path)

    val_confusion_matrix <- confusionMatrix(
      data = val_pred_tbl$pred, reference = val_pred_tbl$obs, mode = "everything", positive = "Responder"
    )
    saveRDS(val_confusion_matrix, val_confusion_matrix_path)
  }
}


# Which cell type to train the model on
which_cell_type <- "CD8 T"
which_cell_type <- "Mono"
which_cell_type <- "all"

# Create train and validation set
partition_by <- "sequencing"
partition_by <- "donor"
if ( partition_by == "sequencing" ) {
  # 1. Partitioned by sequencing methods
  cells_train <- pbmc_int@meta.data %>% dplyr::filter(orig.ident == "PMLCITEseq") %>% rownames()
  cells_val <- pbmc_int@meta.data %>% dplyr::filter(orig.ident == "PMLRNA") %>% rownames()
  saving_flag <- stringr::str_replace(which_cell_type, " ", "_") %>% paste0(".partition_by_", partition_by)
  train_val_loop(pbmc_int, cells_train, cells_val, saving_flag, object_dir)
} else if (partition_by == "donor") {
  # 2. Partitioned by random random donors
  responder_samples <- pbmc_int@meta.data %>% dplyr::select(Response, Vireo_assignment) %>%
    dplyr::distinct() %>% dplyr::filter(Response == "Responder") %>% dplyr::pull(Vireo_assignment)
  non_responder_samples <- pbmc_int@meta.data %>% dplyr::select(Response, Vireo_assignment) %>%
    dplyr::distinct() %>% dplyr::filter(Response == "Non-responder") %>% dplyr::pull(Vireo_assignment)

  set.seed(31415926)
  seed_list <- sample(1000000:100000000, 36)
  seed_idx <- 0
  for (per_responder in responder_samples) {
    for (per_non_responder in non_responder_samples) {
      seed_idx <- seed_idx + 1
      set.seed(seed_list[seed_idx])
      val_donors <- c(per_responder, per_non_responder)
      cells_val <- pbmc_int@meta.data %>% dplyr::filter(Vireo_assignment %in% val_donors) %>% rownames()
      cells_train <- pbmc_int@meta.data %>% dplyr::filter(! Vireo_assignment %in% val_donors) %>% rownames() %>% sample(length(cells_val))
      saving_flag <- stringr::str_replace(which_cell_type, " ", "_") %>% paste0(".partition_by_", partition_by, ".excl_") %>%
        paste0(., paste0(val_donors, collapse = "_"))
      train_val_loop(pbmc_int, cells_train, cells_val, saving_flag, object_dir)
    }
  }

  # Collect metric results
  for (which_metric in c("train", "validation")) {
    all_metrics <- NULL
    for (per_responder in responder_samples) {
      for (per_non_responder in non_responder_samples) {
        val_donors <- c(per_responder, per_non_responder)
        loading_flag <- stringr::str_replace(which_cell_type, " ", "_") %>%
          paste0(".partition_by_", partition_by, ".excl_") %>%
          paste0(., paste0(val_donors, collapse = "_"))
        in_file <- file.path(object_dir, paste0(which_metric, ".", loading_flag, ".gbm_eval.eval_metric.rds"))
        metric <- readRDS(in_file)
        all_metrics <- metric$stdres$Group1 %>%
          tibble::rownames_to_column("metric") %>%
          dplyr::mutate(responder = per_responder, non_responder = per_non_responder) %>%
          dplyr::bind_rows(all_metrics)
      }
    }
    auc_roc_tbl <- all_metrics %>% dplyr::filter(metric == "AUC-ROC") %>%
      dplyr::arrange(Score) %>%
      dplyr::mutate(
        CI_lower = stringr::str_extract(CI, "^(.+)-(.+)$", group = 1),
        CI_upper = stringr::str_extract(CI, "^(.+)-(.+)$", group = 2),
        CI_lower = as.numeric(CI_lower), CI_upper = as.numeric(CI_upper),
        Rank = dplyr::row_number()
      )

    y_min <- auc_roc_tbl$CI_lower %>% min()
    y_max <- auc_roc_tbl$CI_upper %>% max()
    p_rangepoint <- ggplot(auc_roc_tbl) +
      geom_line(aes(x = Rank, y = Score)) +
      geom_ribbon(aes(x = Rank, ymin = CI_lower, ymax = CI_upper), alpha = 0.05) +
      geom_pointrange(aes(x = Rank, y = Score, ymin = CI_lower, ymax = CI_upper, color = Score)) +
      scale_color_gradient(low = "lightblue", high = "darkblue") +
      ylim(y_min, y_max) +
      labs(x = "Rank", y = NULL, color = "AUC-ROC\nScore") +
      theme_classic()

    p_boxplot <- ggplot(data = auc_roc_tbl) +
      geom_violin(aes(x = "AUC-ROC", y = Score), width = 0.75) +
      geom_boxplot(aes(x = "AUC-ROC", y = Score), outlier.shape = NA, width = 0.5) +
      geom_point(aes(x = "AUC-ROC", y = Score), position = "jitter", size = 5, alpha = 0.75, color = "darkblue") +
      ylim(y_min, y_max) +
      labs(x = stringr::str_to_title(which_metric), y = "AUC-ROC") +
      theme_classic()

    p <- (p_boxplot | p_rangepoint) + plot_layout(guides = "collect", widths = c(1, 2))
    plot_save_to <- file.path(plot_dir, paste0(which_metric, ".partition_by_", partition_by, ".eval_metric.pdf"))
    ggsave(plot_save_to, plot = p, width = 7, height = 4)
  }

  # Plotting confusion metric
  selected_params <- c("Accuracy", "Precision", "Recall", "F1", "Kappa", "Prevalence", "Sensitivity", "Specificity")
  for (which_metric in c("train", "validation")) {
    all_metrics <- NULL
    for (per_responder in responder_samples) {
      for (per_non_responder in non_responder_samples) {
        val_donors <- c(per_responder, per_non_responder)
        loading_flag <- stringr::str_replace(which_cell_type, " ", "_") %>%
          paste0(".partition_by_", partition_by, ".excl_") %>%
          paste0(., paste0(val_donors, collapse = "_"))
        in_file <- file.path(object_dir, paste0(which_metric, ".", loading_flag, ".gbm_eval.confusion_matrix.rds"))
        metric <- readRDS(in_file)
        all_metrics <- data.frame(Metric = names(metric$overall), Score = metric$overall) %>%
          dplyr::bind_rows(data.frame(Metric = names(metric$byClass), Score = metric$byClass)) %>%
          dplyr::mutate(responder = per_responder, non_responder = per_non_responder) %>%
          dplyr::bind_rows(all_metrics) %>%
          as.data.table() %>%
          dplyr::filter(Metric %in% selected_params)
      }
    }

    metric_order <- all_metrics %>%
      dplyr::group_by(Metric) %>%
      dplyr::summarise(Score = mean(Score)) %>%
      dplyr::arrange(Score) %>%
      dplyr::pull(Metric)

    p <- ggplot(data = all_metrics) +
      geom_violin(aes(x = Metric, y = Score), width = 0.75) +
      geom_boxplot(aes(x = Metric, y = Score), outlier.shape = NA, width = 0.25) +
      geom_point(aes(x = Metric, y = Score, color = Metric), position = "jitter", size = 1.5, alpha = 0.75) +
      scale_x_discrete(limits = metric_order) +
      labs(x = NULL, y = stringr::str_to_title(which_metric), fill = NULL, color = NULL) +
      theme_classic()
    plot_save_to <- file.path(plot_dir, paste0(which_metric, ".partition_by_", partition_by, ".confusion_matrix.pdf"))
    ggsave(plot_save_to, plot = p, width = 7, height = 4)
  }

  all_metrics <- NULL
  for (per_responder in responder_samples) {
    for (per_non_responder in non_responder_samples) {
      val_donors <- c(per_responder, per_non_responder)
      loading_flag <- stringr::str_replace(which_cell_type, " ", "_") %>%
        paste0(".partition_by_", partition_by, ".excl_") %>%
        paste0(., paste0(val_donors, collapse = "_"))

      in_file <- file.path(object_dir, paste0("train.", loading_flag, ".gbm_eval.variable_importance.rds"))
      metric <- readRDS(in_file)
      all_metrics <- metric$importance %>% as.data.frame() %>%
        tibble::rownames_to_column("Feature") %>%
        dplyr::filter(Overall > 0) %>%
        dplyr::mutate(responder = per_responder, non_responder = per_non_responder) %>%
        dplyr::bind_rows(all_metrics)
    }
  }

  feature_tbl <- all_metrics %>% dplyr::left_join(auc_roc_tbl, by = c("responder", "non_responder"))
  selected_features <- feature_tbl %>% dplyr::group_by(Feature) %>% dplyr::summarise(n = n()) %>% dplyr::filter(n >= 18) %>% dplyr::pull(Feature)

  imp_var_list <- feature_tbl %>%
    dplyr::group_by(responder, non_responder) %>%
    dplyr::summarise(gene_set = list(Feature)) %>%
    dplyr::mutate(gene_set_name = dplyr::row_number()) %>%
    dplyr::pull(gene_set, gene_set_name)
  shared_features <- imp_var_list %>% Reduce(intersect, .)

  plot_tab <- all_metrics %>%
    dplyr::filter(Feature %in% selected_features) %>%
    dplyr::mutate(Model = paste0(responder, ".vs.", non_responder)) %>%
    tidyr::pivot_wider(id_cols = c(Model), names_from = Feature, values_from = Overall) %>%
    dplyr::mutate(dplyr::across(-c(Model), ~ ifelse(is.na(.), 0, .))) %>%
    as.data.frame() %>%
    column_to_rownames("Model") %>%
    as.matrix() %>%
    t()

  rownames(plot_tab)[(rowSums(plot_tab > 10) > 5)]

  marker_features <- c(
    "CCL5", "VIM", "IL32", "IFITM2", "CD74", "MALAT1", "CD52", "NFKB1", "B2M", "FTH1", "HLA-C", "TPT1", "FOS", "EIF1",
    "TMSB10", "FTH1", "MALAT1", "TPT1", "B2M", "EIF1", "DDX5", "S100A4", "NFKB1", "`HLA-A`", "`HLA-C`", "VIM", "BTG1",
    "UBC", "NFKBIA", "FOS", "UBA52", "FTL", "ATP5F1E", "LINC00486", "IL32", "CD247", "`HLA-DRB1`", "IFITM2", "SIPA1L1"
  ) %>% unique
  marker_features <- c(
    "LINC00486", "MALAT1", "DDX5", "NFKB1", "IL32", "IFITM2",  "CD247",
    "CCL5", "VIM", "CD74", "CD52", "B2M"#, "FTH1"
  )

  mark_idx <- match(marker_features, rownames(plot_tab))
  plot_save_to <- file.path(plot_dir, paste0("train.feature_importance.partition_by_", partition_by, ".pdf"))
  pdf(plot_save_to, width = 6, height = 7.5)
  row_ann <- rowAnnotation(
    foo = anno_mark(at = mark_idx, labels = marker_features, labels_gp = gpar(col = "darkred", fontsize = 10, fontface = "italic"))
  )
  Heatmap(
    plot_tab, col = col_fun, name = "Importance",
    right_annotation = row_ann,
    show_row_names = FALSE, show_column_names = FALSE,
  ) %>% draw()
  dev.off()
  # ggsave(plot_save_to, plot = p, width = 7, height = 4)


  ego <- enrichGO(
    gene = rownames(plot_tab), OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "ALL", pAdjustMethod = "BH",
    pvalueCutoff = 0.05, qvalueCutoff = 0.05
  )
  ego_save_to <- file.path(object_dir, paste0("train.feature_importance.partition_by_", partition_by, ".ego.rds"))
  saveRDS(ego, ego_save_to)
  p <- dotplot(ego, showCategory = 15) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  plot_save_to <- file.path(plot_dir, paste0("train.feature_importance.partition_by_", partition_by, ".ego.pdf"))
  ggsave(plot_save_to, plot = p, width = 8.5, height = 5.5)

  ekegg <- bitr(rownames(plot_tab), "SYMBOL", "ENTREZID", OrgDb = org.Hs.eg.db) %>% dplyr::pull(ENTREZID) %>%
    enrichKEGG(gene = ., organism = "hsa", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
  kegg_save_to <- file.path(object_dir, paste0("train.feature_importance.partition_by_", partition_by, ".kegg.rds"))
  saveRDS(ekegg, kegg_save_to)
  p <- dotplot(ekegg, showCategory = 15) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  plot_save_to <- file.path(plot_dir, paste0("train.feature_importance.partition_by_", partition_by, ".kegg.pdf"))
  ggsave(plot_save_to, plot = p, width = 7, height = 6)



  DefaultAssay(pbmc_int) <- "SCT"
  selected_cells <- pbmc_int@meta.data %>%
    dplyr::filter(predicted.celltype.l1 %in% c("Mono", "B", "NK", "CD8 T", "CD4 T"), Timepoint == "BL") %>% rownames()

  p <- lapply(marker_features, function(x) {
    VlnPlot(pbmc_int[, selected_cells], features = x, split.by = "Response", group.by = "predicted.celltype.l1") +
      labs(x = NULL, y = NULL) +
      scale_color_discrete(breaks = c("Responder", "Non-Responder"), labels = c("Rs", "NR")) +
      theme(title = element_text(size = 10))
  }) %>%
    wrap_plots(ncol = 4, guides = "collect") & theme(legend.position = "bottom")

  plot_save_to <- file.path(plot_dir, paste0("train.feature_importance.partition_by_", partition_by, ".violin_plot.pdf"))
  ggsave(plot_save_to, plot = p, width = 8, height = 5.5)
}

# Temporary
pbmc_cite <- readRDS("/home/zzhang/Documents/projects/wp_pml/outputs/analysis/CITE_seq/objects/integrated/pbmc.cite_seq.integrated.pca_umap_clustered.annotated.rds")
pbmc_cite$Response <- pbmc_cite$Response %>% factor(., levels = c("Responder", "Non-responder"))

selected_cells <- pbmc_cite@meta.data %>% as.data.frame() %>% dplyr::filter(Timepoint %in% "BL") %>% rownames()
p <- FeaturePlot(pbmc_cite[, selected_cells], features = c("IFNG"), split.by = "Response", order = TRUE, ncol = 1)
plot_save_to <- file.path("/home/zzhang/Documents/projects/wp_pml/outputs/analysis/cell_cell_communication/plots/ccc.pml_citeseq.baseline.IFNG.vlnplot.pdf")
ggsave(plot_save_to, plot = p, width = 8, height = 4)
