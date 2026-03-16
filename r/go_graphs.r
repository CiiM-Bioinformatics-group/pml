#!/usr/bin/env Rscript
# File: go_graphs.r
# Author: Zhenhua Zhang
# E-mail: zhenhua.zhang217@gmail.com
# Created: Jan 12, 2024
# Updated: Jan 12, 2024


# GO structure
# inverse_of, transitive_over
# is_a, disjoint_from
#
# relationship
#   - ends_during (happens once)
#   - happens_during
#   - has_part
#   - negatively_regulates
#   - occurs_in
#   - part_of
#   - positively_regulates
#   - regulates
#
# intersection_of
#   - happens_during
#   - has_part
#   - negatively_regulates
#   - occurs_in
#   - part_of 
#   - positively_regulates
#   - regulates

library(tidyverse)
library(igraph)
library(tidygraph)
library(ggraph)
library(ggsci)


go_graphml_saveto <- "/home/zzhang/Documents/projects/resources/GeneOntology/go.graphml"
if (!file.exists(go_graphml_saveto)) {
  # Loading all GO terms
  lines <- readLines("/home/zzhang/Documents/projects/resources/GeneOntology/purl.obolibrary.org/obo/go.obo") %>% purrr::discard(~.x == "")
  las_position <- which(lines %in% "[Typedef]") %>% min() %>% `-`(1)

  # A list of GO terms ready to construct the graph
  vertex_list <- data.frame(from = which(lines %in% "[Term]")) %>%
    dplyr::mutate(to = lead(from - 1, default = las_position)) %>%
    apply(1, function(vec) {
      # vec <- c("from" = 32, "to" = 39)
      per_vertex <- lines[vec["from"]:vec["to"]] %>%
        purrr::keep(~ ! .x %in% c("[Term]", "")) %>%
        stringr::str_split(": ", n = 2, simplify = TRUE) %>%
        as.data.frame() %>%
        dplyr::group_by(V1) %>%
        dplyr::summarize(V2 = list(V2)) %>%
        dplyr::pull(V2, V1)

      per_vertex
    }) %>%
    purrr::set_names(lapply(., function(x) x$id)) %>%
    purrr::keep(~is.null(.x$is_obsolete))

  # A blueprint of the GO graph
  go_graph <- make_empty_graph(length(vertex_list))
  vertex_labels <- names(vertex_list)
  V(go_graph)$name <- vertex_labels

  # Adding edges if two GO term vertices are connected due to "is_a", "part_of" or "regulates"
  for (head_idx in seq_along(vertex_list)) {
    per_vertex <- vertex_list[head_idx][[1]]
    go_graph <- set_vertex_attr(go_graph, "go_term", head_idx, per_vertex$name)
    go_graph <- set_vertex_attr(go_graph, "namespace", head_idx, per_vertex$namespace)
    vertex_is_a <- per_vertex$is_a %>% stringr::str_extract_all("GO:[0-9]+") %>% unlist()
    if (!is.null(vertex_is_a)) {
      for (tail_goid in vertex_is_a[vertex_is_a %in% vertex_labels]) {
        tail_idx <- which(vertex_labels == tail_goid)
        go_graph <- add_edges(go_graph, c(head_idx, tail_idx), attr = list("is_a" = TRUE))
      }
    } else {
      cat("[W]: No is_a for", per_vertex$id, "\r")
    }

    vertex_relationship <- per_vertex$relationship
    if (!is.null(vertex_relationship)) {
      cur_relationship <- vertex_relationship %>% stringr::str_split(" ", n = 3)
      for (per_rel in cur_relationship) {
        h2t_rel <- per_rel[1]
        tail_goid <- per_rel[2]
        tail_idx <- which(vertex_labels == tail_goid)

        x2y_edge <- get.edge.ids(go_graph, c(head_idx, tail_idx))
        if (x2y_edge == 0) {
          go_graph <- add_edges(go_graph, c(head_idx, tail_idx), attr = list("relationship" = h2t_rel))
        } else {
          go_graph <- set_edge_attr(go_graph, "relationship", x2y_edge, h2t_rel)
        }
      }
    } else {
      cat("[W]: No relationship for", per_vertex$id, "\r")
    }

    if (head_idx == length(vertex_list)) { cat("\n") }
  }

  write.graph(go_graph, file = go_graphml_saveto, format = "graphml")
} else {
  go_graph <- read.graph(go_graphml_saveto, format = "graphml")
}


go_id <- "GO:0000102"
vs_idx <- as_ids(V(go_graph)) %>% head(500)
sub_go_graph <- subgraph(go_graph, vs_idx)
sub_go_graph_tbl <- as_tbl_graph(sub_go_graph)
sub_go_graph_tbl %>% tidygraph::filter(!is.na(is_a))

g <- ggraph(sub_go_graph_tbl, layout = 'kk') +
  geom_node_point(aes(colour = namespace, size = centrality_pagerank())) +
  geom_edge_link(aes(edge_colour = factor(relationship)), start_cap = circle(0.01, 'npc'), end_cap = circle(0.01, 'npc'), arrow = grid::arrow(length = unit(0.01, 'npc'), type="closed"), check_overlap = TRUE) +
  geom_edge_link(aes(edge_linetype = factor(is_a)), start_cap = circle(0.01, 'npc'), end_cap = circle(0.01, 'npc'), arrow = grid::arrow(length = unit(0.01, 'npc'), type = "open")) +
  scale_edge_linetype(na.value = "blank") +
  scale_edge_colour_discrete(na.value = "white") +
  scale_color_npg() +
  labs(color = "Namespace", edge_colour = "Relationship", edge_linetype = "Is_a", size = "PageRank") +
  theme(legend.position = "right")
ggsave("/vol/projects/BIIM/PML/temp/go-basic.pdf", g, width = 7.5, height = 6)
