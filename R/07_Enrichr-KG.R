# Script name: 07_Enrichr_KG.R
# ==============================================================================
# Title: miRNA target network visualization from Enrichr-KG datasets.

# Author: Sergio Pérez Oliveira

# Description: This script loads and visualizes an interaction network of miRNAs and their target genes, 
#              integrating enrichment results from multiple sources (e.g., GO, KEGG, DisGeNET). It provides 
#              both an overview of the full network and a filtered subnetwork focusing on HTT, ATXN1, and 
#              ATXN2 and their immediate neighbors. Nodes are annotated and colored by source, and the network 
#              is visualized using both static and interactive layouts with `visNetwork`. 
# 
# Inputs: 
#   - nodes.tsv: node metadata including id, label, kind, and color (optional)
#   - edges.tsv: edge list including source, target, and relation
#
# Outputs:
#   - Interactive network plots highlighting key target genes and node types
#   - A subnetwork plot restricted to HTT, ATXN1, ATXN2 and their 1st-degree neighbors
#
# Dependencies:
#   - igraph
#   - tidyverse
#   - ggraph
#   - tidygraph
#   - visNetwork
# ==============================================================================

# Load data ----

nodes <- read.delim(file.choose())   # Select your nodes file
edges <- read.delim(file.choose())   # Select your edges file

head(nodes)
head(edges)

# Initial visualization with visNetwork ----

# Prepare nodes for visNetwork
nodes_vis <- nodes[, c("id", "label", "kind", "color")]
nodes_vis$color <- ifelse(is.na(nodes_vis$color), "gray", nodes_vis$color)  # default to gray if missing

# Prepare edges for visNetwork
edges_vis <- edges
colnames(edges_vis)[colnames(edges_vis) == "source"] <- "from"
colnames(edges_vis)[colnames(edges_vis) == "target"] <- "to"
edges_vis$title <- edges_vis$relation  # shows relation on hover

# Plot interactive network
visNetwork(nodes_vis, edges_vis) %>%
  visPhysics(
    solver = "forceAtlas2Based",
    forceAtlas2Based = list(gravitationalConstant = -150),
    stabilization = list(enabled = TRUE, iterations = 1000),
    minVelocity = 0.1
  ) %>%
  visNodes(   
    color = list(border = "black", background = nodes_vis$color, highlight = "yellow"),
    font = list(size = 16, face = "bold", vadjust = 0),
    scaling = list(min = 10, max = 30)
  ) %>%
  visEdges(smooth = FALSE) %>%
  visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
  visLayout(randomSeed = 123)

# Highlight target genes (HTT, ATXN1, ATXN2) ----

# Step 1: define target genes and color for each
target_genes <- c("HTT", "ATXN1", "ATXN2")
target_colors <- c("HTT" = "red", "ATXN1" = "blue", "ATXN2" = "green")

# Step 2: get IDs of the target nodes
target_node_ids <- nodes_vis$id[nodes_vis$label %in% target_genes]

# Step 3: initialize edge and node styles
edges_vis$color <- "gray"
edges_vis$width <- 1
edges_vis$dashes <- FALSE
nodes_vis$color_current <- nodes_vis$color

# Step 4: update color and width for each target gene
for (gene in target_genes) {
  gene_id <- nodes_vis$id[nodes_vis$label == gene]
  
  # Update node color
  nodes_vis$color_current[nodes_vis$id == gene_id] <- target_colors[gene]
  
  # Update edges connected to this gene
  idx_edges <- which(edges_vis$from == gene_id | edges_vis$to == gene_id)
  edges_vis$color[idx_edges] <- target_colors[gene]
  edges_vis$width[idx_edges] <- 3
}

# Step 5: visualize with updated styles
visNetwork(nodes_vis, edges_vis) %>%
  visPhysics(
    solver = "forceAtlas2Based",
    forceAtlas2Based = list(gravitationalConstant = -150),
    stabilization = TRUE
  ) %>%
  visNodes(
    color = list(border = "black", background = nodes_vis$color_current),
    font = list(size = 16, face = "bold", vadjust = 0, align = "center"),
    scaling = list(min = 10, max = 40)
  ) %>%
  visEdges(smooth = FALSE) %>%
  visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
  visLayout(randomSeed = 123)


# Filter network to target genes and neighbors ----

# Rename columns to match igraph requirements
colnames(edges)[1:2] <- c("from", "to")
colnames(nodes)[1] <- "id"

# Create full graph
graph <- graph_from_data_frame(d = edges, vertices = nodes, directed = TRUE)

# Assign labels
V(graph)$label <- nodes$label[match(V(graph)$name, nodes$id)]

# Select target genes
target_nodes <- V(graph)[V(graph)$label %in% target_genes]

# Get 1st-degree neighbors (ego network)
neighbors <- ego(graph, order = 1, nodes = target_nodes, mode = "all")
included_nodes <- unique(unlist(neighbors))
subgraph <- induced_subgraph(graph, vids = included_nodes)

# Assign labels in the subgraph
V(subgraph)$label <- V(graph)$label[match(V(subgraph)$name, V(graph)$name)]


# Assign colors by node type ----

# Define color palette
node_types <- unique(V(subgraph)$kind)
palette_colors <- c(
  "DisGeNET"                            = "#D2DB7D",
  "GO Biological Process 2021"          = "#FFBAFB",
  "Gene"                                = "#C5E1A5",
  "Human Phenotype Ontology"            = "#B6D7FF",
  "KEGG 2021 Human"                     = "#E9E7F0",
  "MGI Mammalian Phenotype Level 4 2021"= "#FF9600"
)

# Prepare subgraph nodes
sub_nodes <- data.frame(
  id = V(subgraph)$name,
  label = V(subgraph)$label,
  kind = V(subgraph)$kind,
  color = palette_colors[V(subgraph)$kind],
  stringsAsFactors = FALSE
)

# Prepare subgraph edges
sub_edges <- igraph::as_data_frame(subgraph, what = "edges")

# Add legend for node types ----

legend_nodes <- data.frame(
  label = names(palette_colors),
  shape = "dot",
  color = unname(palette_colors),
  size = 15,
  stringsAsFactors = FALSE
)

# Simplify legend labels
legend_nodes$label <- c("DisGeNET", "GO BP", "Gene", "HP", "KEGG", "MP")


# Final visualization of subgraph ----

visNetwork(sub_nodes, sub_edges) %>%
  visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
  visNodes(
    color = list(border = "black", background = sub_nodes$color, highlight = "yellow"),
    borderWidth = 1,
    font = list(face = "bold")
  ) %>%
  visLayout(randomSeed = 123) %>%
  visPhysics(enabled = FALSE) %>%
  visLegend(
    useGroups = FALSE,
    addNodes = legend_nodes,
    position = "left"
  )
