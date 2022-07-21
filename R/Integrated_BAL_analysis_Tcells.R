# R script for scRNA-seq analysis of integrated BAL data Mono/Macro

# seed for reproducibility
seed <- 1999

# Load dataset of merged T-cells
ds <- readRDS("path/to/object")

# Normalize data and identify highly variable features
ds <- Seurat::NormalizeData(
  object               = ds,
  assay                = "RNA",
  normalization.method = "LogNormalize",
  scale.factor         = 10000
)
ds <- Seurat::FindVariableFeatures(
  object           = ds,
  selection.method = "vst", 
  nfeatures        = 3000
)

ds <- Seurat::ScaleData(
  object          = ds
)

ds <- Seurat::RunPCA(
  object = ds,
  npcs   = 55
)

# MNN-corrected PCA, based on Nomralized data 
set.seed(seed = seed)
ds@reductions[["mnn"]] <- Seurat::CreateDimReducObject(
  embeddings = batchelor::fastMNN(
    ds@assays$RNA@data[ds@assays$RNA@var.features, ],
    batch = ds@meta.data$batch,
    d     = 55
  )@int_colData$reducedDims$corrected,
  key        = "MNN_",
  assay      = "RNA"
)

# UMAP with PCA
ds <- Seurat::RunUMAP(
  object         = ds,
  dims           = 1:55,
  seed.use       = seed,
  reduction      = "pca",
  reduction.name = "pca_umap"
)

# UMAP with MNN-corrected PCA; return.model = TRUE for projection later
ds <- Seurat::RunUMAP(
  object         = ds,
  dims           = 1:55,
  seed.use       = seed,
  return.model   = TRUE,
  reduction      = "mnn",
  reduction.name = "mnn_umap")

# Clustering
set.seed(seed = seed)
ds <- Seurat::FindNeighbors(
  object    = ds,
  dims      = 1:55,
  reduction = "mnn"
)

ds <- Seurat::FindClusters(
  object     = ds,
  resolution = 1.5,
  random.seed = seed
)


# Plot Umaps from PCA and MNN-corrected PCA with batch
# Plot PCA_UMAP with  batch
embedding <- "pca_umap"
data <- tidyr::as_tibble(
  x = ds@reductions[[embedding]]@cell.embeddings
)
colnames(data) <- c("x", "y")
category <- "batch"
data$col <- as.factor(ds@meta.data[[category]])
plot <- ggplot2::ggplot(
  data    = data,
  mapping = ggplot2::aes(
    x   = x,
    y   = y,
    col = col
  )
) +
  ggplot2::geom_point(
    size = 0.3
  ) +
  ggplot2::theme(
    panel.background = ggplot2::element_blank(),
    axis.ticks       = ggplot2::element_blank(),
    axis.text        = ggplot2::element_blank(),
    axis.title       = ggplot2::element_text(
      size = 15, hjust = 0.1, vjust = 0.1
    ),
    legend.text      = ggplot2::element_text(size = 15),
    legend.title     = ggplot2::element_blank(),
    legend.position  = "right",
    legend.key       = ggplot2::element_rect(fill = "white", color = "white")
  ) +
  ggplot2::guides(
    color = ggplot2::guide_legend(
      override.aes = list(size = 4),
      reverse = FALSE
    )) +
  ggplot2::coord_fixed() +
  ggplot2::labs(x = "UMAP 1", y = "UMAP 2") +
  ggplot2::annotate(
    geom  = "segment", 
    x     = min(data$x)*1.1, 
    xend  = seq(min(data$x), max(data$x), length.out = 100)[30], 
    y     = min(data$y)*1.1, 
    yend  = min(data$y)*1.1, 
    arrow = ggplot2::arrow(angle = 25, type = "closed")
  ) +
  ggplot2::annotate(
    geom  = "segment", 
    x     = min(data$x)*1.1, 
    xend  = min(data$x)*1.1, 
    y     = min(data$y)*1.1, 
    yend  = seq(min(data$y), max(data$y), length.out = 100)[30], 
    arrow = ggplot2::arrow(angle = 25, type = "closed")
  )
ggplot2::ggsave(
  plot = plot,
  filename = paste0(plot_path, "PCA_Umap_batch.png"),
  width = 9,
  height = 7
)


# Plot MNN-corrected PCA_UMAP with batch
embedding <- "mnn_umap"
data <- tidyr::as_tibble(
  x = ds@reductions[[embedding]]@cell.embeddings
)
colnames(data) <- c("x", "y")
category <- "batch"
data$col <- as.factor(ds@meta.data[[category]])
plot <- ggplot2::ggplot(
  data    = data,
  mapping = ggplot2::aes(
    x   = x,
    y   = y,
    col = col
  )
) +
  ggplot2::geom_point(
    size = 0.3
  ) +
  ggplot2::theme(
    panel.background = ggplot2::element_blank(),
    axis.ticks       = ggplot2::element_blank(),
    axis.text        = ggplot2::element_blank(),
    axis.title       = ggplot2::element_text(
      size = 15, hjust = 0.1, vjust = 0.1
    ),
    legend.text      = ggplot2::element_text(size = 15),
    legend.title     = ggplot2::element_blank(),
    legend.position  = "right",
    legend.key       = ggplot2::element_rect(fill = "white", color = "white")
  ) +
  ggplot2::guides(
    color = ggplot2::guide_legend(
      override.aes = list(size = 4),
      reverse = FALSE
    )) +
  ggplot2::coord_fixed() +
  ggplot2::labs(x = "UMAP 1", y = "UMAP 2") +
  ggplot2::annotate(
    geom  = "segment", 
    x     = min(data$x)*1.1, 
    xend  = seq(min(data$x), max(data$x), length.out = 100)[30], 
    y     = min(data$y)*1.1, 
    yend  = min(data$y)*1.1, 
    arrow = ggplot2::arrow(angle = 25, type = "closed")
  ) +
  ggplot2::annotate(
    geom  = "segment", 
    x     = min(data$x)*1.1, 
    xend  = min(data$x)*1.1, 
    y     = min(data$y)*1.1, 
    yend  = seq(min(data$y), max(data$y), length.out = 100)[30], 
    arrow = ggplot2::arrow(angle = 25, type = "closed")
  )
ggplot2::ggsave(
  plot = plot,
  filename = paste0(plot_path, "MNN_Umap_batch.png"),
  width = 12,
  height = 7
)

# Plot MNN-corrected PCA_UMAP with clusters
embedding <- "mnn_umap"
data <- tidyr::as_tibble(
  x = ds@reductions[[embedding]]@cell.embeddings
)
colnames(data) <- c("x", "y")
category <- "seurat_clusters"
data$col <- as.factor(ds@meta.data[[category]])
plot <- ggplot2::ggplot(
  data    = data,
  mapping = ggplot2::aes(
    x   = x,
    y   = y,
    col = col
  )
) +
  ggplot2::geom_point(
    size = 0.3
  ) +
  ggplot2::theme(
    panel.background = ggplot2::element_blank(),
    axis.ticks       = ggplot2::element_blank(),
    axis.text        = ggplot2::element_blank(),
    axis.title       = ggplot2::element_text(
      size = 15, hjust = 0.1, vjust = 0.1
    ),
    legend.text      = ggplot2::element_text(size = 15),
    legend.title     = ggplot2::element_blank(),
    legend.position  = "right",
    legend.key       = ggplot2::element_rect(fill = "white", color = "white")
  ) +
  ggplot2::guides(
    color = ggplot2::guide_legend(
      override.aes = list(size = 4),
      reverse = FALSE
    )) +
  ggplot2::coord_fixed() +
  ggplot2::labs(x = "UMAP 1", y = "UMAP 2") +
  ggplot2::annotate(
    geom  = "segment", 
    x     = min(data$x)*1.1, 
    xend  = seq(min(data$x), max(data$x), length.out = 100)[30], 
    y     = min(data$y)*1.1, 
    yend  = min(data$y)*1.1, 
    arrow = ggplot2::arrow(angle = 25, type = "closed")
  ) +
  ggplot2::annotate(
    geom  = "segment", 
    x     = min(data$x)*1.1, 
    xend  = min(data$x)*1.1, 
    y     = min(data$y)*1.1, 
    yend  = seq(min(data$y), max(data$y), length.out = 100)[30], 
    arrow = ggplot2::arrow(angle = 25, type = "closed")
  )
ggplot2::ggsave(
  plot = plot,
  filename = paste0(plot_path, "MNN_Umap_Clusters.png"),
  width = 12,
  height = 7
)

# UMAP with gray and CP_BAL_T clusters
embedding <- "mnn_umap"
data <- tidyr::as_tibble(
  x = ds@reductions[[embedding]]@cell.embeddings
)
colnames(data) <- c("x", "y")
category <- "CP_Ann"
data$col <- as.factor(ds@meta.data[[category]])
plot <- ggplot2::ggplot(
  data    = data[data$col != "other", ],
  mapping = ggplot2::aes(
    x   = x,
    y   = y,
    col = col
  )
) +
  ggplot2::geom_point(
    data = data,
    ggplot2::aes(x, y),
    color = "gray60", size = 3) +
  ggplot2::geom_point(
    size = 2.5
  ) +
  ggplot2::theme(
    panel.background = ggplot2::element_blank(),
    axis.ticks       = ggplot2::element_blank(),
    axis.text        = ggplot2::element_blank(),
    axis.title       = ggplot2::element_text(
      size = 15, hjust = 0.1, vjust = 0.1
    ),
    legend.text      = ggplot2::element_text(size = 15),
    legend.title     = ggplot2::element_blank(),
    legend.position  = "right",
    legend.key       = ggplot2::element_rect(fill = "white", color = "white")
  ) +
  ggplot2::guides(
    color = ggplot2::guide_legend(
      override.aes = list(size = 4),
      reverse = FALSE
    )) +
  ggplot2::coord_fixed() +
  ggplot2::labs(x = "UMAP 1", y = "UMAP 2") +
  ggplot2::annotate(
    geom  = "segment", 
    x     = min(data$x)*1.1, 
    xend  = seq(min(data$x), max(data$x), length.out = 100)[30], 
    y     = min(data$y)*1.1, 
    yend  = min(data$y)*1.1, 
    arrow = ggplot2::arrow(angle = 25, type = "closed")
  ) +
  ggplot2::annotate(
    geom  = "segment", 
    x     = min(data$x)*1.1, 
    xend  = min(data$x)*1.1, 
    y     = min(data$y)*1.1, 
    yend  = seq(min(data$y), max(data$y), length.out = 100)[30], 
    arrow = ggplot2::arrow(angle = 25, type = "closed")
  )
ggplot2::ggsave(
  plot = plot,
  filename = paste0(plot_path, "MNN_Umap_CPAnn.png"),
  width = 12,
  height = 7
)


# Th17/1-Scores
gene_ids_Nathan <- to_ids[c("CCR6", "DPP4", "CTSH", "LGALS3", "BHLHE40", "PDE4D")]
gene_ids_Gideon <- to_ids[c("RORA", "RORC", "RBPJ", "BHLHE40", "FURIN", "COTL1", "CCL20", "CCR6",
                            "IL23R", "CXCR3", "IFNG", "TNF")]
gene_idsh1 <- to_ids[c("IFNG", "TBX21", "CXCR3", "NKG7", "FASLG",
                         "IL2", "TNF", "IL12RB1", "IL12RB2", "GZMB")]
gene_idsh17 <- to_ids[c("RORC", "IL17A", "IL17F", "IL22", "CCR6", "IL23R", "CTSH", "BASP1")]


CT <- list()
CT$Score_Nathan <- gene_ids_Nathan
CT$Score_Gideon <- gene_ids_Gideon
CT$Th1Score <- gene_idsh1
CT$Th17Score <- gene_idsh17

ds <- Seurat::AddModuleScore(
  object = ds, features = CT, name = names(CT), nbin = 24
)

# Plot Score_Nathan on UMAP
data <- data.frame("cell" = rownames(ds@meta.data),
                   "x" = ds@reductions$mnn_umap@cell.embeddings[,1],
                   "y" = ds@reductions$mnn_umap@cell.embeddings[,2],
                   "Condition" = ds@meta.data$Condition,
                   "Score" = ds@meta.data$Score_Nathan1)
data <- data[order(data$Score),]
limits <- c(-0.5, 1)
data$Score[data$Score > max(limits)] <- max(limits)
data$Score[data$Score < min(limits)] <- min(limits)
plot <- ggplot2::ggplot(
  data    = data,
  mapping = ggplot2::aes(
    x   = x,
    y   = y,
    col = Score
  )
) +
  ggplot2::geom_point(
    size = 0.5
  ) +
  ggplot2::theme(
    panel.background = ggplot2::element_blank(),
    axis.ticks       = ggplot2::element_blank(),
    axis.text        = ggplot2::element_blank(),
    axis.title       = ggplot2::element_text(size = 15, hjust = 0.1, vjust = 0.1),
    legend.text      = ggplot2::element_text(size = 15),
    legend.title     = ggplot2::element_blank(),
    legend.position  = "right",
    legend.key       = ggplot2::element_rect(fill = "white", color = "white"),
    plot.title            = ggplot2::element_text(face = "italic", size = 17, hjust = 0.5)
  ) +
  ggplot2::scale_color_distiller(
    palette = "RdYlBu", direction = -1, limits = limits
  ) +
ggplot2::guides(
  col = ggplot2::guide_colorbar(
    barwidth = 1, barheight = 15, 
    frame.colour = "black", ticks.colour = "black")
  ) +
  ggplot2::coord_fixed() +
  ggplot2::labs(x = "UMAP 1", y = "UMAP 2") +
  ggplot2::annotate(
    geom  = "segment", 
    x     = min(data$x)*1.1, 
    xend  = seq(min(data$x), max(data$x), length.out = 100)[30], 
    y     = min(data$y)*1.1, 
    yend  = min(data$y)*1.1, 
    arrow = ggplot2::arrow(angle = 25, type = "closed")
  ) +
  ggplot2::annotate(
    geom  = "segment", 
    x     = min(data$x)*1.1, 
    xend  = min(data$x)*1.1, 
    y     = min(data$y)*1.1, 
    yend  = seq(min(data$y), max(data$y), length.out = 100)[30], 
    arrow = ggplot2::arrow(angle = 25, type = "closed")
  )
ggplot2::ggsave(
  plot = plot,
  filename = paste0(plot_path, "Score_Nathan_UMAP.png"),
  width = 7,
  height = 3.5
)

# ViolinPlot Nathan Th17/1 Score ~Clusters
data <- data.frame("Value" = ds@meta.data$Score_Nathan1)
data$col <- data$Value
data$col[data$col > max(limits)] <- max(limits)
data$col[data$col < min(limits)] <- min(limits)
data$cluster <- ds@meta.data$seurat_clusters
limits <- c(-0.5, 1)
# Calculate p-values
ps <- data.frame("cluster" = as.character(levels(data$cluster)))
ps$p.value <- NA
for (i in unique(ps$cluster)) {
  ps$p.value[ps$cluster == i] <- wilcox.test(
    x = data$Value[data$cluster == i],
    y = data$Value,
    alternative = "greater"
  )$p.value
}
ps$p.adjust <- p.adjust(ps$p.value)
ps$label <- round(-log10(ps$p.adjust))
ps$label[ps$label == Inf] <- 150


plot <- ggplot2::ggplot(
  data    = data,
  mapping = ggplot2::aes(
    x     = cluster,
    y     = Value,
    col   = col,
    fill  = cluster
  )
) +
  ggplot2::geom_point(
    position = "jitter", size = 2
  ) +
  ggplot2::scale_color_distiller(
    palette = "RdYlBu", direction = -1, limits = limits
  ) +
  ggplot2::geom_hline(yintercept = mean(data$Value), size = 2, color = "gray55") +
  ggplot2::geom_violin(
    scale = "width", size = 1, alpha = 0.9, draw_quantiles = 0.5, color = "black",
    trim = TRUE
  ) +
  ggplot2::theme_classic(base_size = 15) +
  ggplot2::theme(
    legend.position  = "",
    axis.title       = ggplot2::element_blank(),
    axis.text.x      = ggplot2::element_text(),
    strip.background = ggplot2::element_blank()
  ) +
  ggplot2::scale_y_continuous(limits = c(min(data$Value), max(data$Value + 0.5))) +
  ggplot2::labs(col = NULL, x = NULL) +
  ggnewscale::new_scale("fill") +
  ggplot2::geom_label(
    data = ps,
    mapping = ggplot2::aes(
      label = label, y = max(data$Value + 0.3), fill = label, col = NULL
    ),
    size = 4,
    color = "white"
  ) +
  viridis::scale_fill_viridis(direction = -1, option = "A")

ggplot2::ggsave(
  plot = plot,
  filename = paste0(plot_path, "Score_Nathan_Violin.png"),
  width = 6,
  height = 3.5
)

# Plot Score_Gideon on UMAP
limits <- c(-0.4, 0.6)
data <- data.frame("cell" = rownames(ds@meta.data),
                   "x" = ds@reductions$mnn_umap@cell.embeddings[,1],
                   "y" = ds@reductions$mnn_umap@cell.embeddings[,2],
                   "Condition" = ds@meta.data$Condition,
                   "Score" = ds@meta.data$Score_Gideon2)
data <- data[order(data$Score),]
data$Score[data$Score > max(limits)] <- max(limits)
data$Score[data$Score < min(limits)] <- min(limits)
plot <- ggplot2::ggplot(
  data    = data,
  mapping = ggplot2::aes(
    x   = x,
    y   = y,
    col = Score
  )
) +
  ggplot2::geom_point(
    size = 0.5
  ) +
  ggplot2::theme(
    panel.background = ggplot2::element_blank(),
    axis.ticks       = ggplot2::element_blank(),
    axis.text        = ggplot2::element_blank(),
    axis.title       = ggplot2::element_text(size = 15, hjust = 0.1, vjust = 0.1),
    legend.text      = ggplot2::element_text(size = 15),
    legend.title     = ggplot2::element_blank(),
    legend.position  = "right",
    legend.key       = ggplot2::element_rect(fill = "white", color = "white"),
    plot.title            = ggplot2::element_text(face = "italic", size = 17, hjust = 0.5)
  ) +
  ggplot2::scale_color_distiller(
    palette = "RdYlBu", direction = -1, limits = limits
  ) +
  ggplot2::guides(
    col = ggplot2::guide_colorbar(
      barwidth = 1, barheight = 15, 
      frame.colour = "black", ticks.colour = "black")
  ) +
  ggplot2::coord_fixed() +
  ggplot2::labs(x = "UMAP 1", y = "UMAP 2") +
  ggplot2::annotate(
    geom  = "segment", 
    x     = min(data$x)*1.1, 
    xend  = seq(min(data$x), max(data$x), length.out = 100)[30], 
    y     = min(data$y)*1.1, 
    yend  = min(data$y)*1.1, 
    arrow = ggplot2::arrow(angle = 25, type = "closed")
  ) +
  ggplot2::annotate(
    geom  = "segment", 
    x     = min(data$x)*1.1, 
    xend  = min(data$x)*1.1, 
    y     = min(data$y)*1.1, 
    yend  = seq(min(data$y), max(data$y), length.out = 100)[30], 
    arrow = ggplot2::arrow(angle = 25, type = "closed")
  )
ggplot2::ggsave(
  plot = plot,
  filename = paste0(plot_path, "Score_Gideon_UMAP.png"),
  width = 7,
  height = 3.5
)

# ViolinPlot Gideon Th17/1 Score ~Clusters
limits <- c(-0.4, 0.6)
data <- data.frame("Value" = ds@meta.data$Score_Gideon2)
data$col <- data$Value
data$col[data$col > max(limits)] <- max(limits)
data$col[data$col < min(limits)] <- min(limits)
data$cluster <- ds@meta.data$seurat_clusters

# Calculate p-values
ps <- data.frame("cluster" = as.character(levels(data$cluster)))
ps$p.value <- NA
for (i in unique(ps$cluster)) {
  ps$p.value[ps$cluster == i] <- wilcox.test(
    x = data$Value[data$cluster == i],
    y = data$Value,
    alternative = "greater"
  )$p.value
}
ps$p.adjust <- p.adjust(ps$p.value)
ps$label <- round(-log10(ps$p.adjust))
ps$label[ps$label == Inf] <- 150


plot <- ggplot2::ggplot(
  data    = data,
  mapping = ggplot2::aes(
    x     = cluster,
    y     = Value,
    col   = col,
    fill  = cluster
  )
) +
  ggplot2::geom_point(
    position = "jitter", size = 2
  ) +
  ggplot2::scale_color_distiller(
    palette = "RdYlBu", direction = -1, limits = limits
  ) +
  ggplot2::geom_hline(yintercept = mean(data$Value), size = 2, color = "gray55") +
  ggplot2::geom_violin(
    scale = "width", size = 1, alpha = 0.9, draw_quantiles = 0.5, color = "black",
    trim = TRUE
  ) +
  ggplot2::theme_classic(base_size = 15) +
  ggplot2::theme(
    legend.position  = "",
    axis.title       = ggplot2::element_blank(),
    axis.text.x      = ggplot2::element_text(),
    strip.background = ggplot2::element_blank()
  ) +
  ggplot2::scale_y_continuous(limits = c(min(data$Value), max(data$Value + 0.5))) +
  ggplot2::labs(col = NULL, x = NULL) +
  ggnewscale::new_scale("fill") +
  ggplot2::geom_label(
    data = ps,
    mapping = ggplot2::aes(
      label = label, y = max(data$Value) + 0.3, fill = label, col = NULL
    ),
    size = 4,
    color = "white"
  ) +
  viridis::scale_fill_viridis(direction = -1, option = "A")
ggplot2::ggsave(
  plot = plot,
  filename = paste0(plot_path, "Score_Gideon_Violin.png"),
  width = 8.3,
  height = 5
)


# Gene Heatmap----------------------------
data <- data.frame("Subset" = ds@meta.data$T_subset,
                   "CD4" = ds@assays$RNA@data[to_ids["CD4"],],
                   "CD40LG" = ds@assays$RNA@data[to_ids["CD40LG"],],
                   "CD8A" = ds@assays$RNA@data[to_ids["CD8A"],],
                   "CD8B" = ds@assays$RNA@data[to_ids["CD8B"],],
                   "IFNG" = ds@assays$RNA@data[to_ids["IFNG"],],
                   "TNF" = ds@assays$RNA@data[to_ids["TNF"],],
                   "TBX21" = ds@assays$RNA@data[to_ids["TBX21"],],
                   "CXCR3" = ds@assays$RNA@data[to_ids["CXCR3"],],
                   "RORC" = ds@assays$RNA@data[to_ids["RORC"],],
                   "CCR6" = ds@assays$RNA@data[to_ids["CCR6"],],
                   "KLRB1" = ds@assays$RNA@data[to_ids["KLRB1"],],
                   "CCL20" = ds@assays$RNA@data[to_ids["CCL20"],],
                   "IL23R" = ds@assays$RNA@data[to_ids["IL23R"],],
                   "IL4I1" = ds@assays$RNA@data[to_ids["IL4I1"],],
                   "ABCB1" = ds@assays$RNA@data[to_ids["ABCB1"],],
                   "CSF2" = ds@assays$RNA@data[to_ids["CSF2"],],
                   "CXCR6" = ds@assays$RNA@data[to_ids["CXCR6"],],
                   "IL7R" = ds@assays$RNA@data[to_ids["IL7R"],],
                   "SELL" = ds@assays$RNA@data[to_ids["SELL"],],
                   "CCR4" = ds@assays$RNA@data[to_ids["CCR4"],],
                   "S1PR1" = ds@assays$RNA@data[to_ids["S1PR1"],],
                   "TCF7" = ds@assays$RNA@data[to_ids["TCF7"],],
                   "CCR7" = ds@assays$RNA@data[to_ids["CCR7"],],
                   "FOXP3" = ds@assays$RNA@data[to_ids["FOXP3"],],
                   "IL2RA" = ds@assays$RNA@data[to_ids["IL2RA"],],
                   "SLC4A10" = ds@assays$RNA@data[to_ids["SLC4A10"],],
                   "ME1" = ds@assays$RNA@data[to_ids["ME1"],],
                   "CCL4" = ds@assays$RNA@data[to_ids["CCL4"],],
                   "CCL5" = ds@assays$RNA@data[to_ids["CCL5"],],
                   "NKG7" = ds@assays$RNA@data[to_ids["NKG7"],],
                   "PRF1" = ds@assays$RNA@data[to_ids["PRF1"],],
                   "FCGR3A" = ds@assays$RNA@data[to_ids["FCGR3A"],],
                   "ISG15" = ds@assays$RNA@data[to_ids["ISG15"],],
                   "MX1" = ds@assays$RNA@data[to_ids["MX1"],]
)
library(dplyr)
data <- data %>% group_by(Subset) %>% summarise(across(everything(), list(mean)))
colnames(data) <- c("Subset", gsub("_1", "", colnames(data[, 2:length(colnames(data))])))
data <- as.matrix(data[, 2:length(colnames(data))])
data <- scale(data)
data[data > 1] <- 1
data[data < -1] <- -1

rann <- data.frame("Cluster" = levels(ds@meta.data$T_subset))
rownames(data) <- levels(ds@meta.data$T_subset)
rownames(rann) <- rownames(data)
data <- t(data)
map <- pheatmap::pheatmap(mat = data,
                          color = colorRampPalette(rev(RColorBrewer::brewer.pal(n=7, name = "RdBu")))(500)[1:500],
                          scale = "none",
                          cluster_rows = F,
                          cluster_cols = F,
                          annotation_names_row = TRUE
)
ggplot2::ggsave(
  plot = map,
  filename = paste0(plot_path, "Heatmap.png"),
  width = 3,
  height = 6
)


# Annotate Clusters
cluster.annotation <- c(
  "0"  = "Th17.1",
  "1"  = "Th1-like",
  "2"  = "CD8 CTL",
  "3"  = "CD8 CTL CD16+",
  "4"  = "CD8 CTL CD16+",
  "5"  = "CD4 TCM S1PR1+ CCR4+",
  "6"  = "CD4 CTL",
  "7"  = "Th17.1",
  "8"  = "CD8 CTL CD16+",
  "9"  = "CD4 TCM ribosomal",
  "10" = "Treg",
  "11" = "CD4 IFN-responsive",
  "12" = "Naive-like",
  "13" = "MAIT",
  "14" = "CD8 IFN-responsive"
)

ds@meta.data$T_subset <- factor(
  cluster.annotation[as.character(ds$seurat_clusters)],
  unique(cluster.annotation)[c(1, 6, 2, 5, 7, 10, 9, 8, 11, 3, 4, 12)])


Subset.colors <- c(
  "Th17.1" = "skyblue",
  "CD4 CTL" = "gold3",
  "Th1-like" = "lightgoldenrod",
  "CD4 TCM S1PR1+ CCR4+" = "springgreen2",
  "CD4 TCM ribosomal" = "darkseagreen",
  "Naive-like" = "palegreen",
  "CD4 IFN-responsive" = "plum2",
  "Treg" = "skyblue4",
  "MAIT" = "slateblue",
  "CD8 CTL" = "lightsalmon2",
  "CD8 CTL CD16+" = "lightsalmon4",
  "CD8 IFN-responsive" = "palevioletred2"
)

# Plot MNN-corrected PCA_UMAP with Subsets
embedding <- "mnn_umap"
data <- tidyr::as_tibble(
  x = ds@reductions[[embedding]]@cell.embeddings
)
colnames(data) <- c("x", "y")
category <- "T_subset"
data$col <- as.factor(ds@meta.data[[category]])
plot <- ggplot2::ggplot(
  data    = data,
  mapping = ggplot2::aes(
    x   = x,
    y   = y,
    col = col
  )
) +
  ggplot2::geom_point(
    size = 0.3
  ) +
  ggplot2::scale_color_manual(values = Subset.colors) +
  ggplot2::theme(
    panel.background = ggplot2::element_blank(),
    axis.ticks       = ggplot2::element_blank(),
    axis.text        = ggplot2::element_blank(),
    axis.title       = ggplot2::element_text(
      size = 15, hjust = 0.1, vjust = 0.1
    ),
    legend.text      = ggplot2::element_text(size = 15),
    legend.title     = ggplot2::element_blank(),
    legend.position  = "right",
    legend.key       = ggplot2::element_rect(fill = "white", color = "white")
  ) +
  ggplot2::guides(
    color = ggplot2::guide_legend(
      override.aes = list(size = 4),
      reverse = FALSE
    )) +
  ggplot2::coord_fixed() +
  ggplot2::labs(x = "UMAP 1", y = "UMAP 2") +
  ggplot2::annotate(
    geom  = "segment", 
    x     = min(data$x)*1.1, 
    xend  = seq(min(data$x), max(data$x), length.out = 100)[30], 
    y     = min(data$y)*1.1, 
    yend  = min(data$y)*1.1, 
    arrow = ggplot2::arrow(angle = 25, type = "closed")
  ) +
  ggplot2::annotate(
    geom  = "segment", 
    x     = min(data$x)*1.1, 
    xend  = min(data$x)*1.1, 
    y     = min(data$y)*1.1, 
    yend  = seq(min(data$y), max(data$y), length.out = 100)[30], 
    arrow = ggplot2::arrow(angle = 25, type = "closed")
  )
ggplot2::ggsave(
  plot = plot,
  filename = paste0(plot_path, "MNN_Umap_Subsets.png"),
  width = 12,
  height = 7
)


# ViolinPlot Nathan Th17/1 Score ~Subsets
data <- data.frame("Value" = ds@meta.data$Score_Nathan1)
data$col <- data$Value
data$col[data$col > max(limits)] <- max(limits)
data$col[data$col < min(limits)] <- min(limits)
data$cluster <- ds@meta.data$T_subset
limits <- c(-0.5, 1)
# Calculate p-values
ps <- data.frame("cluster" = as.character(levels(data$cluster)))
ps$p.value <- NA
for (i in unique(ps$cluster)) {
  ps$p.value[ps$cluster == i] <- wilcox.test(
    x = data$Value[data$cluster == i],
    y = data$Value,
    alternative = "greater"
  )$p.value
}
ps$p.adjust <- p.adjust(ps$p.value)
ps$label <- round(-log10(ps$p.adjust))
ps$label[ps$label == Inf] <- 150


plot <- ggplot2::ggplot(
  data    = data,
  mapping = ggplot2::aes(
    x     = cluster,
    y     = Value,
    col   = col,
    fill  = cluster
  )
) +
  ggplot2::geom_point(
    position = "jitter", size = 1.5
  ) +
  ggplot2::scale_color_distiller(
    palette = "RdYlBu", direction = -1, limits = limits
  ) +
  ggplot2::geom_hline(yintercept = mean(data$Value), size = 2, color = "gray55") +
  ggplot2::geom_violin(
    scale = "width", size = 1, alpha = 0.9, draw_quantiles = 0.5, color = "black",
    trim = TRUE
  ) +
  ggplot2::scale_fill_manual(values = Subset.colors) +
  ggplot2::theme_classic(base_size = 15) +
  ggplot2::theme(
    legend.position  = "",
    axis.title       = ggplot2::element_blank(),
    axis.text.x      = ggplot2::element_text(size = 5),
    strip.background = ggplot2::element_blank()
  ) +
  ggplot2::scale_y_continuous(limits = c(min(data$Value), max(data$Value + 0.5))) +
  ggplot2::labs(col = NULL, x = NULL) +
  ggnewscale::new_scale("fill") +
  ggplot2::geom_label(
    data = ps,
    mapping = ggplot2::aes(
      label = label, y = max(data$Value + 0.3), fill = label, col = NULL
    ),
    size = 4,
    color = "white"
  ) +
  viridis::scale_fill_viridis(direction = -1, option = "A")
ggplot2::ggsave(
  plot = plot,
  filename = paste0(plot_path, "Score_Nathan_Violin_Subsets.png"),
  width = 5,
  height = 3
)


# ViolinPlot Gideon Th17/1 Score ~Subsets
data <- data.frame("Value" = ds@meta.data$Score_Gideon2)
data$col <- data$Value
data$col[data$col > max(limits)] <- max(limits)
data$col[data$col < min(limits)] <- min(limits)
data$cluster <- ds@meta.data$T_subset
limits <- c(-0.4, 0.6)
# Calculate p-values
ps <- data.frame("cluster" = as.character(levels(data$cluster)))
ps$p.value <- NA
for (i in unique(ps$cluster)) {
  ps$p.value[ps$cluster == i] <- wilcox.test(
    x = data$Value[data$cluster == i],
    y = data$Value,
    alternative = "greater"
  )$p.value
}
ps$p.adjust <- p.adjust(ps$p.value)
ps$label <- round(-log10(ps$p.adjust))
ps$label[ps$label == Inf] <- 150


plot <- ggplot2::ggplot(
  data    = data,
  mapping = ggplot2::aes(
    x     = cluster,
    y     = Value,
    col   = col,
    fill  = cluster
  )
) +
  ggplot2::geom_point(
    position = "jitter", size = 1.5
  ) +
  ggplot2::scale_color_distiller(
    palette = "RdYlBu", direction = -1, limits = limits
  ) +
  ggplot2::geom_hline(yintercept = mean(data$Value), size = 2, color = "gray55") +
  ggplot2::geom_violin(
    scale = "width", size = 1, alpha = 0.9, draw_quantiles = 0.5, color = "black",
    trim = TRUE
  ) +
  ggplot2::scale_fill_manual(values = Subset.colors) +
  ggplot2::theme_classic(base_size = 15) +
  ggplot2::theme(
    legend.position  = "",
    axis.title       = ggplot2::element_blank(),
    axis.text.x      = ggplot2::element_text(size = 5),
    strip.background = ggplot2::element_blank()
  ) +
  ggplot2::scale_y_continuous(limits = c(min(data$Value), max(data$Value + 0.5))) +
  ggplot2::labs(col = NULL, x = NULL) +
  ggnewscale::new_scale("fill") +
  ggplot2::geom_label(
    data = ps,
    mapping = ggplot2::aes(
      label = label, y = max(data$Value + 0.3), fill = label, col = NULL
    ),
    size = 4,
    color = "white"
  ) +
  viridis::scale_fill_viridis(direction = -1, option = "A")
ggplot2::ggsave(
  plot = plot,
  filename = paste0(plot_path, "Score_Gideon_Violin_Subsets.png"),
  width = 5,
  height = 3
)

# ViolinPlot Gideon Th17/1 Score ~Condition
data <- data.frame("Value" = ds@meta.data$Score_Gideon2)
data$col <- ds@meta.data$Condition
plot <- ggplot2::ggplot(
  data    = data,
  mapping = ggplot2::aes(
    x     = col,
    y     = Value,
    fill  = col
  )
) +
  ggplot2::geom_point(
    position = "jitter", size = 0.3, color = "black", 
  ) +
  ggplot2::geom_hline(yintercept = mean(data$Value), size = 2, color = "gray55") +
  ggplot2::geom_violin(
    scale = "width", size = 1, alpha = 0.75, draw_quantiles = 0.5,
    trim = TRUE
  ) +
  ggplot2::theme_classic(base_size = 15) +
  ggplot2::theme(
    legend.position  = "",
    axis.title       = ggplot2::element_blank(),
    axis.text.x      = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
    strip.background = ggplot2::element_blank()
  )
ggplot2::ggsave(
  plot = plot,
  filename = paste0(plot_path, "Score_Gideon_Violin_Condition.png"),
  width = 5,
  height = 5
)

# ViolinPlot Nathan Th17/1 Score ~Subsets
data <- data.frame("Value" = ds@meta.data$Score_Nathan1)
data$col <- ds@meta.data$Condition
plot <- ggplot2::ggplot(
  data    = data,
  mapping = ggplot2::aes(
    x     = col,
    y     = Value,
    fill  = col
  )
) +
  ggplot2::geom_point(
    position = "jitter", size = 0.3, color = "black", 
  ) +
  ggplot2::geom_hline(yintercept = mean(data$Value), size = 2, color = "gray55") +
  ggplot2::geom_violin(
    scale = "width", size = 1, alpha = 0.75, draw_quantiles = 0.5,
    trim = TRUE
  ) +
  ggplot2::theme_classic(base_size = 15) +
  ggplot2::theme(
    legend.position  = "",
    axis.title       = ggplot2::element_blank(),
    axis.text.x      = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
    strip.background = ggplot2::element_blank()
  )
ggplot2::ggsave(
  plot = plot,
  filename = paste0(plot_path, "Score_Nathan_Violin_Condition.png"),
  width = 5,
  height = 5
)

# Density plot for every Condition
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

data <- data.frame("cell" = rownames(ds@meta.data),
                   "x" = ds@reductions$mnn_umap@cell.embeddings[,1],
                   "y" = ds@reductions$mnn_umap@cell.embeddings[,2],
                   "Condition" = ds@meta.data$Condition)

ggplot2::ggplot(data = data) + ggplot2::geom_point(ggplot2::aes(x, y))

data$density <- NA
data[data$Condition == "COV",]$density <- scales::rescale(get_density(data[data$Condition == "COV",]$x, data[data$Condition == "COV",]$y, n = 100))
data[data$Condition == "CP",]$density <- scales::rescale(get_density(data[data$Condition == "CP",]$x, data[data$Condition == "CP",]$y, n = 100))
data[data$Condition == "HC",]$density <- scales::rescale(get_density(data[data$Condition == "HC",]$x, data[data$Condition == "HC",]$y, n = 100))
data[data$Condition == "Sarc",]$density <- scales::rescale(get_density(data[data$Condition == "Sarc",]$x, data[data$Condition == "Sarc",]$y, n = 100))

plt <- ggplot2::ggplot(data) + ggplot2::geom_point(ggplot2::aes(x, y), color = "gray60", size = 1)


# Condition = Covid19
plot <- ggplot2::ggplot(
  data    = data[data$Condition == "COV", ],
  mapping = ggplot2::aes(
    x   = x,
    y   = y,
    col = density
  )
) +
  ggplot2::geom_point(
    data = data,
    ggplot2::aes(x, y),
    color = "gray60", size = 3) +
  ggplot2::geom_point(
    size = 2
  ) +
  ggplot2::theme(
    panel.background = ggplot2::element_blank(),
    axis.ticks       = ggplot2::element_blank(),
    axis.text        = ggplot2::element_blank(),
    axis.title       = ggplot2::element_text(size = 15, hjust = 0.1, vjust = 0.1),
    legend.text      = ggplot2::element_text(size = 15),
    legend.title     = ggplot2::element_text(face = "italic", size = 18),
    legend.position  = "",
    legend.key       = ggplot2::element_rect(fill = "white", color = "white"),
    plot.title            = ggplot2::element_text(face = "italic", size = 17, hjust = 0.5)
  ) +
  ggplot2::guides(
    color = ggplot2::guide_legend(
      override.aes = list(size = 5),
      reverse = TRUE
    )) +
  ggplot2::scale_color_gradient2(
    low = "white",
    mid = rgb(70, 175, 50, maxColorValue = 255),
    high = rgb(15, 45, 15, maxColorValue = 255),
    midpoint = 0.5) +
  ggplot2::guides(
    color = ggplot2::guide_colorbar(
      barwidth  = 0.6,
      barheight = 10,
      ticks     = TRUE, 
      frame.colour = "black",
      title = "",
      label = FALSE
    )
  ) +
  ggplot2::coord_fixed() +
  ggplot2::labs(x = "UMAP 1", y = "UMAP 2") +
  ggplot2::annotate(
    geom  = "segment", 
    x     = min(data$x)*1.1, 
    xend  = seq(min(data$x), max(data$x), length.out = 100)[30], 
    y     = min(data$y)*1.1, 
    yend  = min(data$y)*1.1, 
    arrow = ggplot2::arrow(angle = 25, type = "closed")
  ) +
  ggplot2::annotate(
    geom  = "segment", 
    x     = min(data$x)*1.1, 
    xend  = min(data$x)*1.1, 
    y     = min(data$y)*1.1, 
    yend  = seq(min(data$y), max(data$y), length.out = 100)[30], 
    arrow = ggplot2::arrow(angle = 25, type = "closed")
  )
ggplot2::ggsave(
  plot = plot,
  filename = paste0(plot_path, "Density_COV.png"),
  width = 8.3,
  height = 5
)

# Condition = CP
plot <- ggplot2::ggplot(
  data    = data[data$Condition == "CP", ],
  mapping = ggplot2::aes(
    x   = x,
    y   = y,
    col = density
  )
) +
  ggplot2::geom_point(
    data = data,
    ggplot2::aes(x, y),
    color = "gray60", size = 3) +
  ggplot2::geom_point(
    size = 2
  ) +
  ggplot2::theme(
    panel.background = ggplot2::element_blank(),
    axis.ticks       = ggplot2::element_blank(),
    axis.text        = ggplot2::element_blank(),
    axis.title       = ggplot2::element_text(size = 15, hjust = 0.1, vjust = 0.1),
    legend.text      = ggplot2::element_text(size = 15),
    legend.title     = ggplot2::element_text(face = "italic", size = 18),
    legend.position  = "",
    legend.key       = ggplot2::element_rect(fill = "white", color = "white"),
    plot.title            = ggplot2::element_text(face = "italic", size = 17, hjust = 0.5)
  ) +
  ggplot2::guides(
    color = ggplot2::guide_legend(
      override.aes = list(size = 5),
      reverse = TRUE
    )) +
  ggplot2::scale_color_gradient2(
    low = "white",
    mid = rgb(70, 175, 50, maxColorValue = 255),
    high = rgb(15, 45, 15, maxColorValue = 255),
    midpoint = 0.5) +
  ggplot2::guides(
    color = ggplot2::guide_colorbar(
      barwidth  = 0.6,
      barheight = 10,
      ticks     = TRUE, 
      frame.colour = "black",
      title = "",
      label = FALSE
    )
  ) +
  ggplot2::coord_fixed() +
  ggplot2::labs(x = "UMAP 1", y = "UMAP 2") +
  ggplot2::annotate(
    geom  = "segment", 
    x     = min(data$x)*1.1, 
    xend  = seq(min(data$x), max(data$x), length.out = 100)[30], 
    y     = min(data$y)*1.1, 
    yend  = min(data$y)*1.1, 
    arrow = ggplot2::arrow(angle = 25, type = "closed")
  ) +
  ggplot2::annotate(
    geom  = "segment", 
    x     = min(data$x)*1.1, 
    xend  = min(data$x)*1.1, 
    y     = min(data$y)*1.1, 
    yend  = seq(min(data$y), max(data$y), length.out = 100)[30], 
    arrow = ggplot2::arrow(angle = 25, type = "closed")
  )
ggplot2::ggsave(
  plot = plot,
  filename = paste0(plot_path, "Density_CP.png"),
  width = 8.3,
  height = 5
)

# Condition = HC
plot <- ggplot2::ggplot(
  data    = data[data$Condition == "HC", ],
  mapping = ggplot2::aes(
    x   = x,
    y   = y,
    col = density
  )
) +
  ggplot2::geom_point(
    data = data,
    ggplot2::aes(x, y),
    color = "gray60", size = 3) +
  ggplot2::geom_point(
    size = 2
  ) +
  ggplot2::theme(
    panel.background = ggplot2::element_blank(),
    axis.ticks       = ggplot2::element_blank(),
    axis.text        = ggplot2::element_blank(),
    axis.title       = ggplot2::element_text(size = 15, hjust = 0.1, vjust = 0.1),
    legend.text      = ggplot2::element_text(size = 15),
    legend.title     = ggplot2::element_text(face = "italic", size = 18),
    legend.position  = "",
    legend.key       = ggplot2::element_rect(fill = "white", color = "white"),
    plot.title            = ggplot2::element_text(face = "italic", size = 17, hjust = 0.5)
  ) +
  ggplot2::guides(
    color = ggplot2::guide_legend(
      override.aes = list(size = 5),
      reverse = TRUE
    )) +
  ggplot2::scale_color_gradient2(
    low = "white",
    mid = rgb(70, 175, 50, maxColorValue = 255),
    high = rgb(15, 45, 15, maxColorValue = 255),
    midpoint = 0.5) +
  ggplot2::guides(
    color = ggplot2::guide_colorbar(
      barwidth  = 0.6,
      barheight = 10,
      ticks     = TRUE, 
      frame.colour = "black",
      title = "",
      label = FALSE
    )
  ) +
  ggplot2::coord_fixed() +
  ggplot2::labs(x = "UMAP 1", y = "UMAP 2") +
  ggplot2::annotate(
    geom  = "segment", 
    x     = min(data$x)*1.1, 
    xend  = seq(min(data$x), max(data$x), length.out = 100)[30], 
    y     = min(data$y)*1.1, 
    yend  = min(data$y)*1.1, 
    arrow = ggplot2::arrow(angle = 25, type = "closed")
  ) +
  ggplot2::annotate(
    geom  = "segment", 
    x     = min(data$x)*1.1, 
    xend  = min(data$x)*1.1, 
    y     = min(data$y)*1.1, 
    yend  = seq(min(data$y), max(data$y), length.out = 100)[30], 
    arrow = ggplot2::arrow(angle = 25, type = "closed")
  )
ggplot2::ggsave(
  plot = plot,
  filename = paste0(plot_path, "Density_HC.png"),
  width = 8.3,
  height = 5
)

# Condition = Sarc
plot <- ggplot2::ggplot(
  data    = data[data$Condition == "Sarc", ],
  mapping = ggplot2::aes(
    x   = x,
    y   = y,
    col = density
  )
) +
  ggplot2::geom_point(
    data = data,
    ggplot2::aes(x, y),
    color = "gray60", size = 3) +
  ggplot2::geom_point(
    size = 2
  ) +
  ggplot2::theme(
    panel.background = ggplot2::element_blank(),
    axis.ticks       = ggplot2::element_blank(),
    axis.text        = ggplot2::element_blank(),
    axis.title       = ggplot2::element_text(size = 15, hjust = 0.1, vjust = 0.1),
    legend.text      = ggplot2::element_text(size = 15),
    legend.title     = ggplot2::element_text(face = "italic", size = 18),
    legend.position  = "",
    legend.key       = ggplot2::element_rect(fill = "white", color = "white"),
    plot.title            = ggplot2::element_text(face = "italic", size = 17, hjust = 0.5)
  ) +
  ggplot2::guides(
    color = ggplot2::guide_legend(
      override.aes = list(size = 5),
      reverse = TRUE
    )) +
  ggplot2::scale_color_gradient2(
    low = "white",
    mid = rgb(70, 175, 50, maxColorValue = 255),
    high = rgb(15, 45, 15, maxColorValue = 255),
    midpoint = 0.5) +
  ggplot2::guides(
    color = ggplot2::guide_colorbar(
      barwidth  = 0.6,
      barheight = 10,
      ticks     = TRUE, 
      frame.colour = "black",
      title = "",
      label = FALSE
    )
  ) +
  ggplot2::coord_fixed() +
  ggplot2::labs(x = "UMAP 1", y = "UMAP 2") +
  ggplot2::annotate(
    geom  = "segment", 
    x     = min(data$x)*1.1, 
    xend  = seq(min(data$x), max(data$x), length.out = 100)[30], 
    y     = min(data$y)*1.1, 
    yend  = min(data$y)*1.1, 
    arrow = ggplot2::arrow(angle = 25, type = "closed")
  ) +
  ggplot2::annotate(
    geom  = "segment", 
    x     = min(data$x)*1.1, 
    xend  = min(data$x)*1.1, 
    y     = min(data$y)*1.1, 
    yend  = seq(min(data$y), max(data$y), length.out = 100)[30], 
    arrow = ggplot2::arrow(angle = 25, type = "closed")
  )
ggplot2::ggsave(
  plot = plot,
  filename = paste0(plot_path, "Density_Sarc.png"),
  width = 8.3,
  height = 5
)


# ViolinPlot Scores in Sarc-Th17.1 cells-------------------------
data <- data.frame(
  "Value" = c(ds_T@meta.data$Score_Nathan1, ds_T@meta.data$Score_Gideon2)
  )
data$Score <- c(rep("1Nathan", ncol(ds_T)), rep("2Gideon", ncol(ds_T)))
data$subset <- c(ds_T@meta.data$T_subset, ds_T@meta.data$T_subset)
data$col <- c(ds_T@meta.data$Condition, ds_T@meta.data$Condition)
# calculate p-values
ps <- data.frame("Score" = c("1Nathan", "2Gideon"))
ps$p.value <- NA
ps$p.value[ps$Score == "1Nathan"] <- wilcox.test(
  x = data$Value[data$Score == "1Nathan" & data$subset == "Th17.1" & data$col == "Sarc"],
  y = data$Value[data$Score == "1Nathan"]
)$p.value
ps$p.value[ps$Score == "2Gideon"] <- wilcox.test(
  x = data$Value[data$Score == "2Gideon" & data$subset == "Th17.1" & data$col == "Sarc"],
  y = data$Value[data$Score == "2Gideon"]
)$p.value

ps$p.adjust <- p.adjust(ps$p.value)
ps$label <- round(-log10(ps$p.adjust))


plot <- ggplot2::ggplot(
  data    = data[data$subset == "Th17.1" & data$col == "Sarc", ],
  mapping = ggplot2::aes(
    x     = Score,
    y     = Value,
    fill  = col
  )
) +
  ggplot2::geom_point(
    position = "jitter", size = 0.3, color = "black", 
  ) +
  ggplot2::geom_segment(ggplot2::aes(
    y = mean(ds_T@meta.data$Score_Gideon2),
    yend = mean(ds_T@meta.data$Score_Gideon2),
    x = 1.5,
    xend = 3),
    size = 2,
    color = "gray55") +
  ggplot2::geom_segment(ggplot2::aes(
    y = mean(ds_T@meta.data$Score_Nathan1),
    yend = mean(ds_T@meta.data$Score_Nathan1),
    x = 0,
    xend = 1.5),
    size = 2,
    color = "gray55") +
  ggplot2::geom_violin(
    scale = "width", size = 1, alpha = 0.75, draw_quantiles = 0.5,
    trim = TRUE
  ) +
  ggplot2::scale_fill_manual(values = c(rgb(red = 207, green = 48, blue = 0, maxColorValue = 255))) +
  ggplot2::theme_classic(base_size = 15) +
  ggplot2::theme(
    legend.position  = "right",
    axis.title       = ggplot2::element_blank(),
    axis.text.x      = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
    strip.background = ggplot2::element_blank()
  ) +
  ggplot2::scale_y_continuous(limits = c(min(data$Value[data$col == "Sarc" & data$subset == "Th17.1"]),
                                         max(data$Value[data$col == "Sarc" & data$subset == "Th17.1"] + 0.5))) +
ggplot2::labs(col = NULL, x = NULL) +
  ggnewscale::new_scale("fill") +
  ggplot2::geom_label(
    data = ps,
    mapping = ggplot2::aes(
      label = label, y = max(data$Value + 0.3), fill = label, col = NULL
    ),
    size = 4,
    color = "white"
  ) +
  viridis::scale_fill_viridis(direction = -1, option = "A")
ggplot2::ggsave(
  plot = plot,
  filename = paste0(plot_path, "Scores_Violin_SarcTh171.png"),
  width = 5,
  height = 5
)
