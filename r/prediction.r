#!/usr/bin/env Rscript
# File: prediction.r
# Author: Zhenhua Zhang
# E-mail: zhenhua.zhang217@gmail.com
# Created: Feb 02, 2024
# Updated: Feb 02, 2024

suppressPackageStartupMessages({
  library(lobstr)

  library(data.table)
  library(tidyverse)
  library(ComplexHeatmap)
  library(clusterProfiler)
  library(ggVennDiagram)
  library(ggsci)
  library(ggbreak)
  library(org.Hs.eg.db)
})

proj_dir <- "~/Documents/projects/wp_pml"

shap_feature_cor_tab_path <- file.path(proj_dir, "outputs/analysis/prediction/function_analysis/correlation_between_shap_and_expression.csv")
shap_feature_cor_tab <- fread(shap_feature_cor_tab_path) %>%
  dplyr::mutate(CellType = dplyr::case_when(
    CellType == "cd4_t" ~ "CD4 T",
    CellType == "cd8_t" ~ "CD8 T",
    CellType == "monocyte" ~ "Mono",
    CellType == "nk_cell" ~ "NK",
    CellType == "b_cell" ~ "B"
)) %>%
  dplyr::mutate(SpearmanFDR < p.adjust(SpearmanPvalue)) %>%
  dplyr::filter(SpearmanFDR < 0.05, abs(SpearmanRho) > 0.1) %>%
  dplyr::select(Feature, CellType, P_val_adj = SpearmanFDR, Statistic = SpearmanRho) %>%
  dplyr::mutate(DataSource = "SHAP")

deg_tab_path <- file.path(proj_dir, "outputs/analysis/deg/pbmc.pml_citeseq.integrated.de_gene.csv")
deg_tab <- fread(deg_tab_path) %>%
  dplyr::filter(donors == "All", celltype %in% c("CD4 T", "CD8 T", "Mono", "B", "NK"), comparison == "Rs.vs.NRs_BL", p_val_adj < 0.05, abs(avg_log2FC) > 0.1) %>%
  dplyr::select(Feature = gene_symbol, CellType = celltype, P_val_adj = p_val_adj, Statistic = avg_log2FC) %>%
  dplyr::mutate(DataSource = "DEG")

# Overlap between DEGs and SHAP-prioritized features
fb_data <- dplyr::left_join(shap_feature_cor_tab, deg_tab, by = c("Feature", "CellType"), suffix = c(".SHAP", ".DEG")) %>%
  dplyr::mutate(FeatureLabel = dplyr::case_when(!is.na(DataSource.DEG) & Statistic.DEG < 0 ~ "Shared(dw)", !is.na(DataSource.DEG) & Statistic.DEG > 0 ~ "Shared(up)", TRUE ~ "Non-shared")) %>%
  dplyr::group_by(CellType, FeatureLabel) %>%
  dplyr::summarise(Count = n()) %>%
  dplyr::group_by(CellType) %>%
  dplyr::mutate(Percent = Count / sum(Count) * 100) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(TotalPercent = Count / sum(Count) * 100) %>%
  dplyr::mutate(Label = paste0(Count, "\n(", round(Percent, 1), "%)")) %>%
  dplyr::group_by(CellType) %>%
  dplyr::mutate(Percent_end = cumsum(Percent), Percent_start = dplyr::lag(Percent_end, default = 0)) %>%
  dplyr::mutate(LabelPos = (Percent_start + Percent_end) / 2) %>%
  dplyr::mutate(LabelPos = dplyr::if_else(LabelPos < 60, 62.5, LabelPos)) %>%
  dplyr::mutate(LabelPos = dplyr::if_else(LabelPos > 97, 100, LabelPos)) %>%
  dplyr::mutate(CellType = factor(CellType, levels = c("CD8 T", "Mono", "CD4 T", "NK", "B"))) %>%
  dplyr::mutate(FeatureLabel = factor(FeatureLabel, levels = c("Shared(up)", "Shared(dw)", "Non-shared")))

fill_plot <- ggplot() +
  geom_col(data = fb_data, aes(x = CellType, y = Percent, fill = FeatureLabel)) +
  geom_text(data = fb_data, aes(x = CellType, y = LabelPos, label = Label), size = 3) +
  scale_fill_jama() +
  scale_y_break(c(0, 60)) +
  labs(x = NULL) +
  theme_classic() +
  theme(legend.position = "top")
fill_plot_path <- file.path(proj_dir, "outputs/analysis/prediction/function_analysis/fill_between_plot.pdf")
ggsave(fill_plot_path, fill_plot, width = 4.75, height = 5)


# Functional enrichment of SHAP features
shap_features <- shap_feature_cor_tab %>%
  dplyr::filter(P_val_adj < 0.05) %>%
  dplyr::group_by(CellType) %>%
  dplyr::summarise(Feature_list = list(Feature)) %>%
  dplyr::pull(Feature_list, CellType)

ego_enr_tab <- lapply(names(shap_features), function(x) {
  ego <- enrichGO(shap_features[[x]], OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "ALL", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
  ego@result %>% as.data.frame() %>% dplyr::mutate(CellType = x)
}) %>%
  dplyr::bind_rows()

ego_plot_tab <- ego_enr_tab %>%
  dplyr::filter(qvalue < 0.05) %>%
  dplyr::group_by(CellType, ONTOLOGY) %>%
  dplyr::slice_min(n = 50, order_by = pvalue) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(CellType = factor(CellType, levels = c("CD8 T", "Mono", "CD4 T", "NK", "B"))) %>%
  dplyr::mutate(Description = forcats::fct_reorder2(Description, Count, CellType))

ego_plot <- ggplot() +
  geom_tile(aes(x = CellType, y = Description, fill = Count), ego_plot_tab) +
  labs(x = "Cell type", y = "GO terms", fill = "Count") +
  theme_classic() +
  theme(legend.position = "top", axis.text.y = element_blank(), axis.ticks.y = element_blank())
ego_plot_path <- file.path(proj_dir, "outputs/analysis/prediction/function_analysis/ego_plot.pdf")
ggsave(ego_plot_path, ego_plot, width = 3, height = 5)


#






# --------- TO BE REMOVED ---------
venn_list = rbind(deg_tab, shap_feature_cor_tab) %>%
  dplyr::summarise(Feature_label = {
    if (n() == 1) {
      return("Non-shared")
    } else {
      if 
    }
  }) %>%
  dplyr::group_by(CellType) %>%
  dplyr::summarise(venn_data = {
    gene_list <- dplyr::cur_data() %>% dplyr::pull(Feature_list, DataSource)
    data.frame(
      DEG = gene_list[["DEG"]][!(gene_list[["DEG"]] %in% gene_list[["SHAP"]])] %>% length(),
      SHAP = gene_list[["SHAP"]][!(gene_list[["SHAP"]] %in% gene_list[["DEG"]])] %>% length(),
      Common = intersect(gene_list[["DEG"]], gene_list[["SHAP"]]) %>% length()
    )}) %>%
  tidyr::unnest(venn_data)
%>%
  dplyr::mutate(List_name = paste0(DataSource, "_", CellType)) %>%
  dplyr::pull(Feature_list, List_name)


# Venn diagram
matrix(names(venn_list), nrow = 5) %>% as.data.frame() %>% apply(1, function(vec) {
  set_a <- venn_list[[vec[1]]]
  set_b <- venn_list[[vec[2]]]
  cell_type <- stringr::str_split(vec[1], "_")[[1]][2]

  set_list <- list(set_a, set_b) %>% purrr::set_names(c("DE\ngenes", "SHAP\nfeatures"))
  venn <- Venn(set_list)
  data <- process_data(venn)
  p <- ggplot() +
    geom_path(aes(X, Y, group = id), color = "black", data = venn_setedge(data), show.legend = FALSE) +
    geom_polygon(aes(X, Y, group = id), fill = "white", alpha = 0.01, data = venn_regionedge(data)) +
    geom_label(aes(X, Y, label = count), data = venn_regionlabel(data)) +
    geom_label(aes(X, Y, label = name), data = venn_setlabel(data)) +
    labs(title = cell_type) +
    coord_fixed(clip = "off") +
    theme_void()

  save_token <- paste0(vec[1], "_vs_", vec[2])
  save_to <- file.path(proj_dir, "outputs/analysis/prediction/function_analysis", paste0("venn_diagram.", save_token, ".pdf"))
  ggsave(save_to, p, width = 2.5, height = 2.75)

  invisible(save_to)
})
