# Script name: 07_Enrichr-KG.R
# ==============================================================================
# Title: Contextual Enrichr-KG network around HTT, ATXN1 and ATXN2.
#
# Author: Sergio Pérez Oliveira
#
# Reviewer-driven interpretation:
#   This analysis is EXPLORATORY and CONTEXTUAL only.
#
#   The Enrichr-KG network is gene-centered and is NOT derived from:
#     - CAG repeat length,
#     - intermediate-allele carrier status,
#     - STR-specific molecular measurements,
#     - repeat-QTL data,
#     - experimentally validated consequences of the repeat variants studied.
#
#   Therefore, network proximity or shared annotations must NOT be interpreted as
#   functional validation, mediation, causality, or a CAG-specific mechanism for
#   the associations observed in HTT, ATXN1 or ATXN2.
#
#   The purpose of this script is limited to showing public-knowledge annotations
#   and first-degree relationships around the three genes as hypothesis-generating
#   biological context.
#
# Inputs:
#   - nodes.tsv: Enrichr-KG node metadata; required columns: id, label, kind.
#                A color column is optional.
#   - edges.tsv: Enrichr-KG edge list; required columns: source, target, relation.
#
# Outputs:
#   - Interactive full Enrichr-KG network.
#   - Interactive 1st-degree contextual subnetwork around HTT/ATXN1/ATXN2.
#   - Machine-readable target-neighbour nodes/edges and network summary in
#     results/reviewer_revision/.
#
# Important:
#   No inferential p-values are generated here. This script must not be used to
#   claim enrichment of CAG-repeat effects or repeat-specific functional evidence.
# ==============================================================================

# Load data ---------------------------------------------------------------------

nodes <- read.delim(
  file.choose(),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

edges <- read.delim(
  file.choose(),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

revision_dir <- file.path("results", "reviewer_revision")
if (!dir.exists(revision_dir)) {
  dir.create(revision_dir, recursive = TRUE)
}

target_genes <- c("HTT", "ATXN1", "ATXN2")

# ==============================================================================
# 1. INPUT VALIDATION AND PROVENANCE SUMMARY
# ==============================================================================

required_node_cols <- c("id", "label", "kind")
required_edge_cols <- c("source", "target", "relation")

missing_node_cols <- setdiff(required_node_cols, names(nodes))
missing_edge_cols <- setdiff(required_edge_cols, names(edges))

if (length(missing_node_cols) > 0) {
  stop(
    "nodes.tsv is missing required column(s): ",
    paste(missing_node_cols, collapse = ", "),
    call. = FALSE
  )
}

if (length(missing_edge_cols) > 0) {
  stop(
    "edges.tsv is missing required column(s): ",
    paste(missing_edge_cols, collapse = ", "),
    call. = FALSE
  )
}

nodes$id <- as.character(nodes$id)
nodes$label <- as.character(nodes$label)
nodes$kind <- as.character(nodes$kind)

edges$source <- as.character(edges$source)
edges$target <- as.character(edges$target)
edges$relation <- as.character(edges$relation)

# Add a default display color if the source file does not provide one.
if (!"color" %in% names(nodes)) {
  nodes$color <- "gray"
}
nodes$color[is.na(nodes$color) | nodes$color == ""] <- "gray"

# Unique node IDs are required for graph construction.
if (anyDuplicated(nodes$id) > 0) {
  duplicate_ids <- unique(nodes$id[duplicated(nodes$id)])
  stop(
    "Duplicate node IDs detected: ",
    paste(utils::head(duplicate_ids, 10), collapse = ", "),
    call. = FALSE
  )
}

# All edge endpoints should refer to known nodes.
known_ids <- nodes$id
unknown_edge_ids <- setdiff(
  unique(c(edges$source, edges$target)),
  known_ids
)

if (length(unknown_edge_ids) > 0) {
  stop(
    "edges.tsv contains endpoint ID(s) absent from nodes.tsv: ",
    paste(utils::head(unknown_edge_ids, 10), collapse = ", "),
    call. = FALSE
  )
}

target_presence <- tibble::tibble(
  gene = target_genes,
  present = target_genes %in% nodes$label
)

if (!all(target_presence$present)) {
  warning(
    "Target gene(s) missing from Enrichr-KG nodes: ",
    paste(target_presence$gene[!target_presence$present], collapse = ", "),
    call. = FALSE
  )
}

node_kind_summary <- nodes %>%
  dplyr::count(kind, name = "n_nodes") %>%
  dplyr::arrange(dplyr::desc(n_nodes))

relation_summary <- edges %>%
  dplyr::count(relation, name = "n_edges") %>%
  dplyr::arrange(dplyr::desc(n_edges))

network_summary <- tibble::tibble(
  metric = c(
    "total_nodes",
    "total_edges",
    "target_genes_requested",
    "target_genes_present"
  ),
  value = c(
    nrow(nodes),
    nrow(edges),
    length(target_genes),
    sum(target_presence$present)
  )
)

readr::write_csv(
  target_presence,
  file.path(revision_dir, "07_target_gene_presence.csv")
)

readr::write_csv(
  network_summary,
  file.path(revision_dir, "07_network_summary.csv")
)

readr::write_csv(
  node_kind_summary,
  file.path(revision_dir, "07_node_kind_summary.csv")
)

readr::write_csv(
  relation_summary,
  file.path(revision_dir, "07_edge_relation_summary.csv")
)

cat("\n============================================================\n")
cat("ENRICHR-KG CONTEXTUAL NETWORK QC\n")
cat("============================================================\n")
print(network_summary, n = Inf, width = Inf)
print(target_presence, n = Inf, width = Inf)

cat("\nInterpretation constraint:\n")
cat(
  "This network is gene-centered public-knowledge context and is not ",
  "CAG/STR-specific functional evidence.\n",
  sep = ""
)

# ==============================================================================
# 2. FULL CONTEXTUAL NETWORK
# ==============================================================================

nodes_vis <- nodes[, unique(c("id", "label", "kind", "color")), drop = FALSE]

edges_vis <- edges
names(edges_vis)[names(edges_vis) == "source"] <- "from"
names(edges_vis)[names(edges_vis) == "target"] <- "to"
edges_vis$title <- edges_vis$relation

# Highlight the three genes without changing the biological meaning of node kind.
nodes_vis$is_target_gene <- nodes_vis$label %in% target_genes
nodes_vis$borderWidth <- ifelse(nodes_vis$is_target_gene, 4, 1)
nodes_vis$size <- ifelse(nodes_vis$is_target_gene, 30, 15)

network_full <- visNetwork::visNetwork(nodes_vis, edges_vis) %>%
  visNetwork::visPhysics(
    solver = "forceAtlas2Based",
    forceAtlas2Based = list(gravitationalConstant = -150),
    stabilization = list(enabled = TRUE, iterations = 1000),
    minVelocity = 0.1
  ) %>%
  visNetwork::visNodes(
    font = list(size = 16, face = "bold", vadjust = 0)
  ) %>%
  visNetwork::visEdges(smooth = FALSE) %>%
  visNetwork::visOptions(
    highlightNearest = TRUE,
    nodesIdSelection = TRUE
  ) %>%
  visNetwork::visLayout(randomSeed = 123)

print(network_full)

# ==============================================================================
# 3. FIRST-DEGREE SUBNETWORK AROUND HTT, ATXN1 AND ATXN2
# ==============================================================================

graph_edges <- edges %>%
  dplyr::transmute(
    from = source,
    to = target,
    relation = relation
  )

graph_nodes <- nodes %>%
  dplyr::rename(name = id)

graph <- igraph::graph_from_data_frame(
  d = graph_edges,
  vertices = graph_nodes,
  directed = TRUE
)

igraph::V(graph)$label <- graph_nodes$label[
  match(igraph::V(graph)$name, graph_nodes$name)
]

target_nodes <- igraph::V(graph)[
  igraph::V(graph)$label %in% target_genes
]

if (length(target_nodes) == 0L) {
  stop(
    "None of HTT, ATXN1 or ATXN2 is present in the graph.",
    call. = FALSE
  )
}

# First-degree neighborhood, irrespective of edge direction, because this is a
# contextual annotation view rather than a causal/directional model.
neighbor_sets <- igraph::ego(
  graph,
  order = 1,
  nodes = target_nodes,
  mode = "all"
)

included_nodes <- unique(unlist(neighbor_sets))
subgraph <- igraph::induced_subgraph(
  graph,
  vids = included_nodes
)

# ==============================================================================
# 4. EXPORT THE EXACT CONTEXT SHOWN
# ==============================================================================

sub_nodes_raw <- igraph::as_data_frame(
  subgraph,
  what = "vertices"
)

sub_edges_raw <- igraph::as_data_frame(
  subgraph,
  what = "edges"
)

sub_nodes_export <- sub_nodes_raw %>%
  dplyr::mutate(
    is_target_gene = label %in% target_genes
  ) %>%
  dplyr::arrange(
    dplyr::desc(is_target_gene),
    kind,
    label
  )

sub_edges_export <- sub_edges_raw %>%
  dplyr::mutate(
    touches_target_gene =
      from %in% sub_nodes_export$name[sub_nodes_export$is_target_gene] |
      to %in% sub_nodes_export$name[sub_nodes_export$is_target_gene]
  )

readr::write_csv(
  sub_nodes_export,
  file.path(revision_dir, "07_target_context_nodes.csv")
)

readr::write_csv(
  sub_edges_export,
  file.path(revision_dir, "07_target_context_edges.csv")
)

context_summary <- tibble::tibble(
  metric = c(
    "context_nodes",
    "context_edges",
    "target_nodes_present",
    "first_degree_non_target_nodes"
  ),
  value = c(
    nrow(sub_nodes_export),
    nrow(sub_edges_export),
    sum(sub_nodes_export$is_target_gene),
    sum(!sub_nodes_export$is_target_gene)
  )
)

readr::write_csv(
  context_summary,
  file.path(revision_dir, "07_target_context_summary.csv")
)

cat("\n============================================================\n")
cat("TARGET-CENTERED CONTEXT SUMMARY\n")
cat("============================================================\n")
print(context_summary, n = Inf, width = Inf)

# ==============================================================================
# 5. TARGET-CENTERED INTERACTIVE VISUALIZATION
# ==============================================================================

palette_colors <- c(
  "DisGeNET"                             = "#D2DB7D",
  "GO Biological Process 2021"           = "#FFBAFB",
  "Gene"                                 = "#C5E1A5",
  "Human Phenotype Ontology"             = "#B6D7FF",
  "KEGG 2021 Human"                      = "#E9E7F0",
  "MGI Mammalian Phenotype Level 4 2021" = "#FF9600"
)

sub_nodes <- sub_nodes_export %>%
  dplyr::transmute(
    id = name,
    label = label,
    kind = kind,
    color = dplyr::coalesce(
      unname(palette_colors[kind]),
      "gray"
    ),
    borderWidth = ifelse(is_target_gene, 4, 1),
    size = ifelse(is_target_gene, 32, 16),
    title = ifelse(
      is_target_gene,
      paste0(
        label,
        "<br><b>Study target gene</b>",
        "<br>Contextual Enrichr-KG annotation only; not CAG-specific evidence."
      ),
      paste0(
        label,
        "<br>Node type: ",
        kind
      )
    )
  )

sub_edges <- sub_edges_export %>%
  dplyr::transmute(
    from = from,
    to = to,
    title = relation
  )

legend_nodes <- data.frame(
  label = c("DisGeNET", "GO BP", "Gene", "HP", "KEGG", "MP"),
  shape = "dot",
  color = unname(palette_colors),
  size = 15,
  stringsAsFactors = FALSE
)

network_target_context <- visNetwork::visNetwork(
  sub_nodes,
  sub_edges
) %>%
  visNetwork::visOptions(
    highlightNearest = TRUE,
    nodesIdSelection = TRUE
  ) %>%
  visNetwork::visNodes(
    font = list(face = "bold")
  ) %>%
  visNetwork::visLayout(
    randomSeed = 123
  ) %>%
  visNetwork::visPhysics(
    enabled = FALSE
  ) %>%
  visNetwork::visLegend(
    useGroups = FALSE,
    addNodes = legend_nodes,
    position = "left"
  )

print(network_target_context)

# ==============================================================================
# 6. MACHINE-READABLE INTERPRETATION NOTE
# ==============================================================================

interpretation_note <- c(
  "07 Enrichr-KG network interpretation",
  "",
  "Status: exploratory / hypothesis-generating only.",
  "",
  "This network is based on gene-centered public knowledge annotations and",
  "first-degree relationships around HTT, ATXN1 and ATXN2.",
  "",
  "It does not use CAG repeat length, intermediate-allele status, STR-specific",
  "functional measurements, repeat-QTL data, or experimental validation of",
  "the repeat variants analyzed in this study.",
  "",
  "Accordingly, network proximity and shared annotations must not be described",
  "as evidence of a CAG-specific molecular mechanism, mediation, causality,",
  "or functional validation of the genetic association results."
)

writeLines(
  interpretation_note,
  con = file.path(
    revision_dir,
    "07_network_interpretation_note.txt"
  )
)

cat("\n============================================================\n")
cat("07 COMPLETE\n")
cat("============================================================\n")
cat(
  "Network retained as contextual/hypothesis-generating analysis only.\n",
  "No CAG-specific functional inference is performed.\n",
  sep = ""
)

# Session info ------------------------------------------------------------------
sessionInfo()
