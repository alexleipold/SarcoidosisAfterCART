
seed <- 1999

ds_all <- readRDS("/home/alexander/Data/members/Alexander_Leipold/projects/Leo_Rasche/MM_Case_Report/Revision/13_05_22/merged_BALs_Annotated.Rds")

ds_M <- subset(ds_all, subset = Celltype %in% c("Mono/Macro", "Mono/Macro prolif"))



# Normalize data and identify highly variable features
ds_M <- Seurat::NormalizeData(
  object               = ds_M,
  assay                = "RNA",
  normalization.method = "LogNormalize",
  scale.factor         = 10000
)
ds_M <- Seurat::FindVariableFeatures(
  object           = ds_M,
  selection.method = "vst", 
  nfeatures        = 1000
)

ds_M <- Seurat::ScaleData(
  object          = ds_M,
  vars.to.regress = "nCount_RNA"
)

ds_M <- Seurat::RunPCA(
  object = ds_M,
  npcs   = 55
)

# MNN-corrected PCA, based on Nomralized data 
set.seed(seed = seed)
ds_M@reductions[["mnn"]] <- Seurat::CreateDimReducObject(
  embeddings = batchelor::fastMNN(
    ds_M@assays$RNA@scale.data[ds_M@assays$RNA@var.features, ],
    batch = ds_M@meta.data$batch,
    d     = 55
  )@int_colData$reducedDims$corrected,
  key        = "MNN_",
  assay      = "RNA"
)

# UMAP with PCA
ds_M <- Seurat::RunUMAP(
  object         = ds_M,
  dims           = 1:25,
  seed.use       = seed,
  reduction      = "pca",
  reduction.name = "pca_umap"
)

# UMAP with MNN-corrected PCA
ds_M <- Seurat::RunUMAP(
  object         = ds_M,
  dims           = 1:25,
  seed.use       = seed,
  reduction      = "mnn",
  reduction.name = "mnn_umap")

# Clustering
set.seed(seed = seed)
ds_M <- Seurat::FindNeighbors(
  object    = ds_M,
  dims      = 1:25,
  reduction = "mnn"
)

ds_M <- Seurat::FindClusters(
  object     = ds_M,
  resolution = 0.4,
  random.seed = seed
)


# Plot PCA_UMAP with batch
embedding <- "pca_umap"
data <- tidyr::as_tibble(
  x = ds_M@reductions[[embedding]]@cell.embeddings
)
colnames(data) <- c("x", "y")
category <- "batch"
data$col <- as.factor(ds_M@meta.data[[category]])
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
  filename = paste0(plot_path, "BAL_Integrated_Mp/pca_Umap_batch.png"),
  width = 8,
  height = 7
)


# Plot MNN-corrected PCA_UMAP with batch
embedding <- "mnn_umap"
data <- tidyr::as_tibble(
  x = ds_M@reductions[[embedding]]@cell.embeddings
)
colnames(data) <- c("x", "y")
category <- "batch"
data$col <- as.factor(ds_M@meta.data[[category]])
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
  filename = paste0(plot_path, "BAL_Integrated_Mp/mnn_Umap_batch.png"),
  width = 8,
  height = 7
)

# Seurat clusters DE analysis
markers_M <- Seurat::FindAllMarkers(ds_M, only.pos = TRUE)
markers_M$gene <- to_genes[markers_M$gene]
View(markers_M)


# Annotate Clusters
cluster.annotation <- c(
  "0"  = "AMp",
  "1"  = "AMp",
  "2"  = "FCN1 Mono",
  "3"  = "CD163/LGMN Mp",
  "4"  = "AMp",
  "5"  = "Mono/Mp",
  "6"  = "Foam cells",
  "7"  = "prolif AMp",
  "8"  = "Pro-inflammatory AMp",
  "9"  = "metallothionein AMp",
  "10" = "LowQ",
  "11" = "IFN-responsive AMp"
)

ds_M@meta.data$Mp_subset <- factor(
  cluster.annotation[as.character(ds_M$seurat_clusters)])

ds_M@meta.data$Mp_subset <- factor(
  cluster.annotation[as.character(ds_M$seurat_clusters)],
  unique(cluster.annotation)[c(2, 4, 3, 1, 7, 8, 10, 5, 6, 9)])



# Density plot for every Condition
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

data <- data.frame("cell" = rownames(ds_M@meta.data),
                   "x" = ds_M@reductions$mnn_umap@cell.embeddings[,1],
                   "y" = ds_M@reductions$mnn_umap@cell.embeddings[,2],
                   "Condition" = ds_M@meta.data$Condition)

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
  filename = paste0(plot_path, "BAL_Integrated_Mp/Density_COV.png"),
  width = 5,
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
  filename = paste0(plot_path, "BAL_Integrated_Mp/Density_CP.png"),
  width = 5,
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
  filename = paste0(plot_path, "BAL_Integrated_Mp/Density_HC.png"),
  width = 5,
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
  filename = paste0(plot_path, "BAL_Integrated_Mp/Density_Sarc.png"),
  width = 5,
  height = 5
)


# Plot Annotation proportions
td <- table(ds_M@meta.data[ds_M@meta.data$Mp_subset != "LowQ",]$batch, ds_M@meta.data[ds_M@meta.data$Mp_subset != "LowQ",]$Mp_subset)
td_to_pct <- function(td) {
  for (i in 1:nrow(td)) {
    pcts <- c()
    for (ii in 1:ncol(td)) {
      pcts <- c(pcts, td[i,ii]/sum(td[i,]))
    }
    td[i, ] <- pcts*100
  }
  return(td)
}
tdpct <- td_to_pct(td)
data <- data.frame(tdpct)
plot <- plot <- ggplot2::ggplot(
  data    = data,
  mapping = ggplot2::aes(
    x     = Var1,
    y     = Freq,
    col   = Var2,
    fill  = Var2
  )
) +
  ggplot2::geom_col(position = "stack") +
  ggplot2::theme(
    panel.background = ggplot2::element_blank(),
    panel.border     = ggplot2::element_rect(size = 2, fill = NA),
    axis.text        = ggplot2::element_text(size = 15),
    axis.text.x      = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1), 
    axis.title       = ggplot2::element_blank(),
    legend.position  = "right",
    legend.text      = ggplot2::element_text(size = 20),
    legend.title     = ggplot2::element_blank()
  ) +
  ggplot2::guides(
    color = ggplot2::guide_legend(
      override.aes = list(label = "")
    )
  )
ggplot2::ggsave(
  plot = plot,
  filename = paste0(plot_path, "BAL_Integrated_Mp/Proportions.png"),
  width = 10,
  height = 6
)


# Plot MNN-corrected PCA_UMAP with Cluster
embedding <- "mnn_umap"
data <- tidyr::as_tibble(
  x = ds_M@reductions[[embedding]]@cell.embeddings
)
colnames(data) <- c("x", "y")
category <- "seurat_clusters"
data$col <- as.factor(ds_M@meta.data[[category]])
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
  filename = paste0(plot_path, "BAL_Integrated_Mp/mnn_Umap_Cluster.png"),
  width = 8,
  height = 7
)

# Plot MNN-corrected PCA_UMAP with Annotation
embedding <- "mnn_umap"
data <- tidyr::as_tibble(
  x = ds_M@reductions[[embedding]]@cell.embeddings
)
colnames(data) <- c("x", "y")
category <- "Mp_subset"
data$col <- as.factor(ds_M@meta.data[[category]])
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
  filename = paste0(plot_path, "BAL_Integrated_Mp/mnn_Umap_Annotation.png"),
  width = 8,
  height = 7
)

# Plot markers as Dotplot
# Dotplot of marker genes
# marker genes
marker_genes <- c(
  "FCN1", "VCAN", "CCR2", # FCN1-Mono
  # Mono-Mp
  "CCL2", "SPP1", "LGMN", # CD163/LGMN Mp
  "FABP4", "INHBA", "CES1", # AMp
  "CCL4", "CCL3", "CCL20", "TNFAIP6", "SOD2", "CXCL8", # Pro inflamma
  "MT1G", "MT2A", "MT1X", "MT1E", "MT1F", # metallo
  "IFIT1", "ISG15", "RSAD2", "MX1", # IFN res
  "LDLR", "SQLE", "HMGCR", "MSMO1", # Foam
  "MKI67", "TOP2A", "CENPA", # prolif
  "CD3E", "IL32", "CD2" # Low Q
)

genes <- rev(marker_genes)
data <- t(as.matrix(ds_M@assays$RNA@data[to_ids[genes], ]))
data <- tidyr::as_tibble(data)
colnames(data) <- genes
data$Cluster <- factor(ds_M@meta.data$seurat_clusters)
data$Cluster <- factor(data$Cluster, levels(data$Cluster)[c(3, 4, 6, 5, 1, 2, 9, 10, 12, 7, 8, 11)])
data$Cluster <- factor(data$Cluster, levels(data$Cluster)[seq(12,1)])
data <- tidyr::gather(data, "Gene", "Expression", -Cluster)
data$Gene <- factor(
  x      = data$Gene,
  levels = rev(unique(data$Gene))
)
data$Pct <- data$Expression > 0
data$Freq <- rep(1, length(data$Pct))
data <- dplyr::group_by(data, Gene)
data <- dplyr::mutate(data, Scaled = scale(Expression)[, 1])
data <- dplyr::group_by(data, Cluster, Gene)
data <- dplyr::summarise(
  data,
  Mean   = mean(Expression),
  Scaled = mean(Scaled),
  Pct    = sum(Pct)/sum(Freq)*100
)
colorscale <- RColorBrewer::brewer.pal(n = 9, name = "RdBu")
data$Scaled[data$Scaled > 1.5] <- 1.5
plot <- ggplot2::ggplot(
  data = data,
  mapping = ggplot2::aes(
    y = Cluster,
    x = Gene,
    size = Pct,
    col = Scaled
  )
) +
  ggplot2::geom_point() +
  ggplot2::scale_color_gradient2(midpoint=0, low=colorscale[length(colorscale)], mid="white",
                                 high=colorscale[1]) +
  ggplot2::theme_classic(base_size = 15) +
  ggplot2::theme(
    axis.title       = ggplot2::element_blank(),
    axis.text.x      = ggplot2::element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      face  = "italic"
    ),
    axis.text.y      = ggplot2::element_text(
      face  = "plain"
    ),
    legend.position  = "right",
    legend.title     = ggplot2::element_blank()
  ) +
  ggplot2::guides(
    color = ggplot2::guide_colorbar(
      barheight = 6, barwidth = 0.6, frame.colour = "black", ticks = FALSE, order = 1
    ),
    size  = ggplot2::guide_legend(
      label.position = "right"
    )
  ) +
  ggplot2::scale_size_area(max_size = 6)
ggplot2::ggsave(
  filename = paste0(plot_path, "BAL_Integrated_Mp/marker_Dotplot.png"),
  width = 8,
  height = 4.5
)




Seurat::DimPlot(ds_M, reduction = "mnn_umap", group.by = "seurat_clusters", split.by = "Condition")

Seurat::FeaturePlot(ds_M, reduction = "mnn_umap", features = to_ids["FABP4"], split.by = "Condition")

Seurat::VlnPlot(ds_M, features = to_ids["IFITM3"])


td <- table(ds_M@meta.data[ds_M@meta.data$Mp_subset != "LowQ",]$batch, ds_M@meta.data[ds_M@meta.data$Mp_subset != "LowQ",]$Mp_subset)

td_to_pct <- function(td) {
  for (i in 1:nrow(td)) {
    pcts <- c()
    for (ii in 1:ncol(td)) {
      pcts <- c(pcts, td[i,ii]/sum(td[i,]))
    }
    td[i, ] <- pcts*100
  }
  return(td)
}
tdpct <- td_to_pct(td)

data <- data.frame(tdpct)

ggplot2::ggplot(data=data, ggplot2::aes(x=Var1, y=Freq, fill=Var2)) +
  ggplot2::geom_bar(stat="identity")

ds_M <- Seurat::SetIdent(ds_M, value = ds_M@meta.data$Mp_subset)
mrkrs_M <- Seurat::FindMarkers(subset(ds_M, subset = Mp_subset == "AMp"), ident.1 = "CP", ident.2 = "Sarc")
mrkrs_M$gene <- to_genes[rownames(mrkrs_M)]
View(mrkrs_M)


