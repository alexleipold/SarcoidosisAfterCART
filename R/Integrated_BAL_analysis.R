# R script for scRNA-seq analysis of integrated BAL data

# set seed for reproducibility
seed <- 1999

# Path for plots
plot_path <- "path/for/plots/"

# Normalize mRNA counts, Identify variable fetures and scale
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

# Calculate low-dimensional embedding & clustering
# PCA
ds <- Seurat::RunPCA(
  object = ds,
  npcs   = 55
  )

# data integration with patient/donor as confounding factor (batch)
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

# UMAP with MNN-corrected PCA
ds <- Seurat::RunUMAP(
  object         = ds,
  dims           = 1:55,
  seed.use       = seed,
  reduction      = "mnn",
  reduction.name = "mnn_umap")

# Plot UMAPS with PCA and MNN-corrected PCA with batch
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
    size = 0.01
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
  width = 10.3,
  height = 9
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
    size = 0.01
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
  height = 9
)


# Downstream MNN-corrected PCA will be used!

# Clustering
set.seed(seed = seed)
ds <- Seurat::FindNeighbors(
  object    = ds,
  dims      = 1:55,
  reduction = "mnn"
  )

ds <- Seurat::FindClusters(
  object     = ds,
  resolution = 0.4
)

# Annotate clusters
cluster.annotation <- c(
  "0"  = "Neutrophil",
  "1"  = "Mono/Macro",
  "2"  = "Mono/Macro",
  "3"  = "T cell",
  "4"  = "Mono/Macro",
  "5"  = "Mono/Macro",
  "6"  = "T cell",
  "7"  = "Neutrophil",
  "8"  = "T cell",
  "9"  = "T prolif",
  "10" = "T prolif",
  "11" = "Plasma",
  "12" = "B/DC",
  "13" = "NK",
  "14" = "Mono/Macro prolif",
  "15" = "Mono/Macro",
  "16" = "Mono/Macro",
  "17" = "Mast"
)

ds@meta.data$Celltype <- factor(
  cluster.annotation[as.character(ds$RNA_snn_res.0.4)],
  unique(cluster.annotation)[c(5, 6, 3, 4, 7, 2, 8, 9, 1)])
                               
# Plot Clusters and Celltypes on UMAP
# Plot Clusters
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
    size = 0.1
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
  filename = paste0(plot_path, "Clusters.png"),
  width = 8.7,
  height = 9
)

# Plot Celltypes
Celltype.colors <- c(
  "T cell" = "skyblue",
  "Mono/Macro" = "mediumorchid2",
  "Mast" = "turquoise",
  "B/DC" = "plum2",
  "NK" = "slateblue",
  "Mono/Macro prolif" = "mediumorchid4",
  "T prolif" = "cyan4",
  "Neutrophil" = "yellow4",
  "Plasma" = "indianred2"
)

embedding <- "mnn_umap"
data <- tidyr::as_tibble(
  x = ds@reductions[[embedding]]@cell.embeddings
)
colnames(data) <- c("x", "y")
category <- "Celltype"
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
    size = 0.1
  ) +
  ggplot2::scale_color_manual(values = Celltype.colors[levels(ds@meta.data$Celltype)]) +
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
  filename = paste0(plot_path, "Celltypes.png"),
  width = 10,
  height = 9
)

# Dotplot of marker genes
# marker genes
marker_genes <- c("JCHAIN", "SDC1", "MZB1", "TNFRSF17", # Plasma
                  "CLEC10A", "CD1C", "MS4A1", "SPIB", # DC/B
                  "CD3D", "CD3E", "IL32", "CD2", # T
                  "MKI67", "TOP2A", # Proliferating
                  "XCL2", "TRDC", "KLRF1", # NK
                  "CD68", "CD14", "CD163", # Mono/Macro
                  "TPSAB1", "MS4A2", "GATA2", # Mast
                  "CXCL8", "CSF3R", "G0S2" # Neutro
                  )

genes <- rev(marker_genes)
data <- t(as.matrix(ds@assays$RNA@data[to_ids[genes], ]))
data <- tidyr::as_tibble(data)
colnames(data) <- genes
data$Cluster <- factor(ds@meta.data$seurat_clusters)
data$Cluster <- factor(data$Cluster, levels(data$Cluster)[c(12, #Plasme
                                                            13, #DC/B
                                                            4, 7, 9, # T
                                                            10, 11, # T prolif
                                                            14, # NK
                                                            2, 3, 5, 6, 16, 17, # Mono/Macro
                                                            15, #Mono/Macro prolif
                                                            18, #Mast
                                                            1, 8 # Neutro
                                                            )])
data$Cluster <- factor(data$Cluster, levels(data$Cluster)[seq(18,1)])
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
data$Scaled[data$Scaled > 2] <- 2
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
  filename = paste0(plot_path, "Dotplot.png"),
  width = 8,
  height = 6
)

saveRDS(ds, file = "/home/alexander/Data/members/Alexander_Leipold/projects/Leo_Rasche/MM_Case_Report/Revision/13_05_22/merged_BALs_Annotated.Rds", compress = FALSE)
