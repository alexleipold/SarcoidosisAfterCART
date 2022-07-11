# R script for scRNA-seq analysis of integrated BAL data Mono/Macro

seed <- 1999


# create gene symbol-EnsemblID lookup table
to_ids <- ds@misc$features$ENSEMBL
names(to_ids) <- ds@misc$features$SYMBOL

to_genes <- ds@misc$features$SYMBOL
names(to_genes) <- ds@misc$features$ENSEMBL

ds <- readRDS("/home/alexander/Data/members/Alexander_Leipold/projects/Leo_Rasche/MM_Case_Report/Revision/13_05_22/merged_BALs_Annotated.Rds")

ds_T <- subset(ds, subset = Celltype %in% c("T cell"))

# get batch in
index <- c("CART_BAL", "CART_BAL",
           "Liao2020_HC1", "Liao2020_HC3",
           "Wendisch2021_B16", "Wendisch2021_B14", "Wendisch2021_B23", "Wendisch2021_B31",
           "Wendisch2021_B20", "Wendisch2021_B29", "Wendisch2021_B28",
           "Liao2021_Sarc1", "Liao2021_Sarc2", "Liao2021_Sarc3", "Liao2021_Sarc4")
names(index) <- unique(ds_T@meta.data$Origin)
ds_T@meta.data$batch <- factor(index[ds_T@meta.data$Origin])
ds_T@meta.data$batch <- factor(ds_T@meta.data$batch,
                             levels(ds_T@meta.data$batch)[c(1, 4, 5, 6, 7, 2, 3, 8, 10, 12, 9, 11, 13, 14)])

# Normalize data and identify highly variable features
ds_T <- Seurat::NormalizeData(
  object               = ds_T,
  assay                = "RNA",
  normalization.method = "LogNormalize",
  scale.factor         = 10000
)
ds_T <- Seurat::FindVariableFeatures(
  object           = ds_T,
  selection.method = "vst", 
  nfeatures        = 3000
)

ds_T <- Seurat::ScaleData(
  object          = ds_T
)

ds_T <- Seurat::RunPCA(
  object = ds_T,
  npcs   = 55
)

# MNN-corrected PCA, based on Nomralized data 
set.seed(seed = seed)
ds_T@reductions[["mnn"]] <- Seurat::CreateDimReducObject(
  embeddings = batchelor::fastMNN(
    ds_T@assays$RNA@data[ds_T@assays$RNA@var.features, ],
    batch = ds_T@meta.data$batch,
    d     = 55
  )@int_colData$reducedDims$corrected,
  key        = "MNN_",
  assay      = "RNA"
)

# UMAP with PCA
ds_T <- Seurat::RunUMAP(
  object         = ds_T,
  dims           = 1:55,
  seed.use       = seed,
  reduction      = "pca",
  reduction.name = "pca_umap"
)

# UMAP with MNN-corrected PCA; return.model = TRUE for projection later
ds_T <- Seurat::RunUMAP(
  object         = ds_T,
  dims           = 1:55,
  seed.use       = seed,
  return.model   = TRUE,
  reduction      = "mnn",
  reduction.name = "mnn_umap")

# Clustering
set.seed(seed = seed)
ds_T <- Seurat::FindNeighbors(
  object    = ds_T,
  dims      = 1:55,
  reduction = "mnn"
)

ds_T <- Seurat::FindClusters(
  object     = ds_T,
  resolution = 1.5,
  random.seed = seed
)


# Plot Umaps from PCA and MNN-corrected PCA with batch
# Plot PCA_UMAP with  batch
embedding <- "pca_umap"
data <- tidyr::as_tibble(
  x = ds_T@reductions[[embedding]]@cell.embeddings
)
colnames(data) <- c("x", "y")
category <- "batch"
data$col <- as.factor(ds_T@meta.data[[category]])
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
  filename = paste0(plot_path, "BAL_Integrated_T/PCA_Umap_batch.png"),
  width = 9,
  height = 7
)


# Plot MNN-corrected PCA_UMAP with batch
embedding <- "mnn_umap"
data <- tidyr::as_tibble(
  x = ds_T@reductions[[embedding]]@cell.embeddings
)
colnames(data) <- c("x", "y")
category <- "batch"
data$col <- as.factor(ds_T@meta.data[[category]])
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
  filename = paste0(plot_path, "BAL_Integrated_T/MNN_Umap_batch.png"),
  width = 12,
  height = 7
)

# Plot MNN-corrected PCA_UMAP with clusters
embedding <- "mnn_umap"
data <- tidyr::as_tibble(
  x = ds_T@reductions[[embedding]]@cell.embeddings
)
colnames(data) <- c("x", "y")
category <- "seurat_clusters"
data$col <- as.factor(ds_T@meta.data[[category]])
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
  filename = paste0(plot_path, "BAL_Integrated_T/MNN_Umap_Clusters.png"),
  width = 12,
  height = 7
)

# UMAP with gray and CP_BAL_T clusters
embedding <- "mnn_umap"
data <- tidyr::as_tibble(
  x = ds_T@reductions[[embedding]]@cell.embeddings
)
colnames(data) <- c("x", "y")
category <- "CP_Ann"
data$col <- as.factor(ds_T@meta.data[[category]])
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
  filename = paste0(plot_path, "BAL_Integrated_T/MNN_Umap_CPAnn.png"),
  width = 12,
  height = 7
)


# Th17/1-Scores
gene_ids_Nathan <- to_ids[c("CCR6", "DPP4", "CTSH", "LGALS3", "BHLHE40", "PDE4D")]
gene_ids_Gideon <- to_ids[c("RORA", "RORC", "RBPJ", "BHLHE40", "FURIN", "COTL1", "CCL20", "CCR6",
                            "IL23R", "CXCR3", "IFNG", "TNF")]
gene_ids_Th1 <- to_ids[c("IFNG", "TBX21", "CXCR3", "NKG7", "FASLG",
                         "IL2", "TNF", "IL12RB1", "IL12RB2", "GZMB")]
gene_ids_Th17 <- to_ids[c("RORC", "IL17A", "IL17F", "IL22", "CCR6", "IL23R", "CTSH", "BASP1")]


CT <- list()
CT$exTh17Score_Nathan <- gene_ids_Nathan
CT$exTh17Score_Gideon <- gene_ids_Gideon
CT$Th1Score <- gene_ids_Th1
CT$Th17Score <- gene_ids_Th17

ds_T <- Seurat::AddModuleScore(
  object = ds_T, features = CT, name = names(CT), nbin = 24
)

# Plot Score_Nathan on UMAP
data <- data.frame("cell" = rownames(ds_T@meta.data),
                   "x" = ds_T@reductions$mnn_umap@cell.embeddings[,1],
                   "y" = ds_T@reductions$mnn_umap@cell.embeddings[,2],
                   "Condition" = ds_T@meta.data$Condition,
                   "Score" = ds_T@meta.data$exTh17Score_Nathan1)
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
  filename = paste0(plot_path, "BAL_Integrated_T/Score_Nathan_UMAP.png"),
  width = 7,
  height = 3.5
)

# ViolinPlot Nathan Th17/1 Score ~Clusters
data <- data.frame("Value" = ds_T@meta.data$exTh17Score_Nathan1)
data$col <- data$Value
data$col[data$col > max(limits)] <- max(limits)
data$col[data$col < min(limits)] <- min(limits)
data$cluster <- ds_T@meta.data$seurat_clusters
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
  filename = paste0(plot_path, "BAL_Integrated_T/Score_Nathan_Violin.png"),
  width = 6,
  height = 3.5
)

# Plot Score_Gideon on UMAP
limits <- c(-0.4, 0.6)
data <- data.frame("cell" = rownames(ds_T@meta.data),
                   "x" = ds_T@reductions$mnn_umap@cell.embeddings[,1],
                   "y" = ds_T@reductions$mnn_umap@cell.embeddings[,2],
                   "Condition" = ds_T@meta.data$Condition,
                   "Score" = ds_T@meta.data$exTh17Score_Gideon2)
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
  filename = paste0(plot_path, "BAL_Integrated_T/Score_Gideon_UMAP.png"),
  width = 7,
  height = 3.5
)

# ViolinPlot Gideon Th17/1 Score ~Clusters
limits <- c(-0.4, 0.6)
data <- data.frame("Value" = ds_T@meta.data$exTh17Score_Gideon2)
data$col <- data$Value
data$col[data$col > max(limits)] <- max(limits)
data$col[data$col < min(limits)] <- min(limits)
data$cluster <- ds_T@meta.data$seurat_clusters

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
  filename = paste0(plot_path, "BAL_Integrated_T/Score_Gideon_Violin.png"),
  width = 8.3,
  height = 5
)


# Gene Heatmap----------------------------
data <- data.frame("Subset" = ds_T@meta.data$T_subset,
                   "CD4" = ds_T@assays$RNA@data[to_ids["CD4"],],
                   "CD40LG" = ds_T@assays$RNA@data[to_ids["CD40LG"],],
                   "CD8A" = ds_T@assays$RNA@data[to_ids["CD8A"],],
                   "CD8B" = ds_T@assays$RNA@data[to_ids["CD8B"],],
                   "IFNG" = ds_T@assays$RNA@data[to_ids["IFNG"],],
                   "TNF" = ds_T@assays$RNA@data[to_ids["TNF"],],
                   "TBX21" = ds_T@assays$RNA@data[to_ids["TBX21"],],
                   "CXCR3" = ds_T@assays$RNA@data[to_ids["CXCR3"],],
                   "RORC" = ds_T@assays$RNA@data[to_ids["RORC"],],
                   "CCR6" = ds_T@assays$RNA@data[to_ids["CCR6"],],
                   "KLRB1" = ds_T@assays$RNA@data[to_ids["KLRB1"],],
                   "CCL20" = ds_T@assays$RNA@data[to_ids["CCL20"],],
                   "IL23R" = ds_T@assays$RNA@data[to_ids["IL23R"],],
                   "IL4I1" = ds_T@assays$RNA@data[to_ids["IL4I1"],],
                   "ABCB1" = ds_T@assays$RNA@data[to_ids["ABCB1"],],
                   "CSF2" = ds_T@assays$RNA@data[to_ids["CSF2"],],
                   "CXCR6" = ds_T@assays$RNA@data[to_ids["CXCR6"],],
                   "IL7R" = ds_T@assays$RNA@data[to_ids["IL7R"],],
                   "SELL" = ds_T@assays$RNA@data[to_ids["SELL"],],
                   "CCR4" = ds_T@assays$RNA@data[to_ids["CCR4"],],
                   "S1PR1" = ds_T@assays$RNA@data[to_ids["S1PR1"],],
                   "TCF7" = ds_T@assays$RNA@data[to_ids["TCF7"],],
                   "CCR7" = ds_T@assays$RNA@data[to_ids["CCR7"],],
                   "FOXP3" = ds_T@assays$RNA@data[to_ids["FOXP3"],],
                   "IL2RA" = ds_T@assays$RNA@data[to_ids["IL2RA"],],
                   "SLC4A10" = ds_T@assays$RNA@data[to_ids["SLC4A10"],],
                   "ME1" = ds_T@assays$RNA@data[to_ids["ME1"],],
                   "CCL4" = ds_T@assays$RNA@data[to_ids["CCL4"],],
                   "CCL5" = ds_T@assays$RNA@data[to_ids["CCL5"],],
                   "NKG7" = ds_T@assays$RNA@data[to_ids["NKG7"],],
                   "PRF1" = ds_T@assays$RNA@data[to_ids["PRF1"],],
                   "FCGR3A" = ds_T@assays$RNA@data[to_ids["FCGR3A"],],
                   "ISG15" = ds_T@assays$RNA@data[to_ids["ISG15"],],
                   "MX1" = ds_T@assays$RNA@data[to_ids["MX1"],]
)
library(dplyr)
data <- data %>% group_by(Subset) %>% summarise(across(everything(), list(mean)))
colnames(data) <- c("Subset", gsub("_1", "", colnames(data[, 2:length(colnames(data))])))
data <- as.matrix(data[, 2:length(colnames(data))])
data <- scale(data)
data[data > 1] <- 1
data[data < -1] <- -1

rann <- data.frame("Cluster" = levels(ds_T@meta.data$T_subset))
rownames(data) <- levels(ds_T@meta.data$T_subset)
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
  filename = paste0(plot_path, "BAL_Integrated_T/Heatmap.png"),
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

ds_T@meta.data$T_subset <- factor(
  cluster.annotation[as.character(ds_T$seurat_clusters)],
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
  x = ds_T@reductions[[embedding]]@cell.embeddings
)
colnames(data) <- c("x", "y")
category <- "T_subset"
data$col <- as.factor(ds_T@meta.data[[category]])
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
  filename = paste0(plot_path, "BAL_Integrated_T/MNN_Umap_Subsets.png"),
  width = 12,
  height = 7
)


# ViolinPlot Nathan Th17/1 Score ~Subsets
data <- data.frame("Value" = ds_T@meta.data$exTh17Score_Nathan1)
data$col <- data$Value
data$col[data$col > max(limits)] <- max(limits)
data$col[data$col < min(limits)] <- min(limits)
data$cluster <- ds_T@meta.data$T_subset
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
  filename = paste0(plot_path, "BAL_Integrated_T/Score_Nathan_Violin_Subsets.png"),
  width = 5,
  height = 3
)


# ViolinPlot Gideon Th17/1 Score ~Subsets
data <- data.frame("Value" = ds_T@meta.data$exTh17Score_Gideon2)
data$col <- data$Value
data$col[data$col > max(limits)] <- max(limits)
data$col[data$col < min(limits)] <- min(limits)
data$cluster <- ds_T@meta.data$T_subset
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
  filename = paste0(plot_path, "BAL_Integrated_T/Score_Gideon_Violin_Subsets.png"),
  width = 5,
  height = 3
)

# ViolinPlot Gideon Th17/1 Score ~Condition
data <- data.frame("Value" = ds_T@meta.data$exTh17Score_Gideon2)
data$col <- ds_T@meta.data$Condition
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
  filename = paste0(plot_path, "BAL_Integrated_T/Score_Gideon_Violin_Condition.png"),
  width = 5,
  height = 5
)

# ViolinPlot Nathan Th17/1 Score ~Subsets
data <- data.frame("Value" = ds_T@meta.data$exTh17Score_Nathan1)
data$col <- ds_T@meta.data$Condition
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
  filename = paste0(plot_path, "BAL_Integrated_T/Score_Nathan_Violin_Condition.png"),
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

data <- data.frame("cell" = rownames(ds_T@meta.data),
                   "x" = ds_T@reductions$mnn_umap@cell.embeddings[,1],
                   "y" = ds_T@reductions$mnn_umap@cell.embeddings[,2],
                   "Condition" = ds_T@meta.data$Condition)

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
  filename = paste0(plot_path, "BAL_Integrated_T/Density_COV.png"),
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
  filename = paste0(plot_path, "BAL_Integrated_T/Density_CP.png"),
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
  filename = paste0(plot_path, "BAL_Integrated_T/Density_HC.png"),
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
  filename = paste0(plot_path, "BAL_Integrated_T/Density_Sarc.png"),
  width = 8.3,
  height = 5
)

# Plot UMAP with genes
pltlist <- list()
genes <- c(
  "IFNG", "TNF", "RORC", "TBX21", "CCR6", "CXCR3", "KLRB1",
  "IL23R", "IL4I1", "ABCB1", "CXCR6", "CSF2", "CCL20", "IL17A",
  
  "CD4", "CD40LG", "CD8A", "CD8B", "STAT3", "RUNX1", "RUNX3", "IL17B", "IL25", "IL17F",
  "IL26", "IL2", "CCR7", "TCF7", "S1PR1", "NKG7", "GZMK", "PRF1", "FOXP3", "IL2RA"
)

for (i in genes) {
  embedding <- "mnn_umap"
  data <- tidyr::as_tibble(
    x = ds_T@reductions[[embedding]]@cell.embeddings
  )
  colnames(data) <- c("x", "y")
  category <- i
  data[["value"]] <- ds_T@assays$RNA@data[to_ids[category],]
  data <- data[order(data[["value"]]), ]
  plot <- ggplot2::ggplot(
    data    = data,
    mapping = ggplot2::aes(
      x   = x,
      y   = y,
      col = value
    )
  ) +
    ggplot2::geom_point(
      size = 0.2
    ) +
    ggplot2::theme(
      panel.background = ggplot2::element_blank(),
      axis.ticks       = ggplot2::element_blank(),
      axis.text        = ggplot2::element_blank(),
      axis.title       = ggplot2::element_blank(),
      legend.text      = ggplot2::element_text(size = 15),
      legend.title     = ggplot2::element_text(face = "italic", size = 18),
      legend.position  = "",
      legend.key       = ggplot2::element_rect(fill = "white", color = "white"),
      plot.title            = ggplot2::element_text(face = "italic", size = 17, hjust = 0.5)
    ) +
    ggplot2::ggtitle(i) +
    ggplot2::guides(
      color = ggplot2::guide_legend(
        override.aes = list(size = 5),
        reverse = TRUE
      )) +
    viridis::scale_color_viridis(
      option    = "B",
      direction = -1
    ) +
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
      arrow = ggplot2::arrow(angle = 25, type = "closed", length = ggplot2::unit(0.1, "inches"))
    ) +
    ggplot2::annotate(
      geom  = "segment", 
      x     = min(data$x)*1.1, 
      xend  = min(data$x)*1.1, 
      y     = min(data$y)*1.1, 
      yend  = seq(min(data$y), max(data$y), length.out = 100)[30], 
      arrow = ggplot2::arrow(angle = 25, type = "closed", length = ggplot2::unit(0.1, "inches"))
    )
  
  assign(paste0("plot_", i), plot)
  pltlist[[paste0("plot_", i)]] <- plot
}


# UMAPxgeneA plots
voidplt <- ggplot2::ggplot() + ggplot2::theme_void()

pltlistFig2 <- list(plot_IFNG, plot_TNF, plot_CSF2, plot_RORC, plot_TBX21, plot_CCR6, plot_CXCR3,
                    voidplt, voidplt, voidplt, voidplt, voidplt, voidplt, voidplt,
                    plot_KLRB1, plot_IL23R, plot_IL4I1, plot_ABCB1, plot_CXCR6, plot_CCL20, plot_IL17A)
plot <- ggpubr::ggarrange(plotlist = pltlistFig2,
                          nrow = 3, ncol = 7,
                          heights = c(1, -0.6, 1))
ggplot2::ggsave(
  plot = plot,
  filename = paste0(plot_path, "BAL_Integrated_T/UMAP_GenesA.png"),
  width = 16,
  height = 8.5
)

# UMAPxgeneB plots
pltlistFigS3 <- list(plot_CD4, plot_CD40LG, plot_CD8A, plot_CD8B, plot_STAT3, plot_RUNX1, plot_RUNX3, plot_IL17B, plot_IL25, plot_IL17F,
                     voidplt, voidplt, voidplt, voidplt, voidplt, voidplt, voidplt, voidplt, voidplt, voidplt,
                     plot_IL26, plot_IL2, plot_CCR7, plot_TCF7, plot_S1PR1, plot_NKG7, plot_GZMK, plot_PRF1, plot_FOXP3, plot_IL2RA)
plot <- ggpubr::ggarrange(plotlist = pltlistFigS3,
                          nrow = 3, ncol = 10,
                          heights = c(1, -0.6, 1))

ggplot2::ggsave(
  plot = plot,
  filename = paste0(plot_path, "BAL_Integrated_T/UMAP_GenesB.png"),
  width = 23,
  height = 8.5
)



### Unifrac---------------
dataall <- as.matrix(cbind(
  subset(ds_T, Condition == "CP")@assays$RNA@data,
  subset(ds_T, Condition == "HC")@assays$RNA@data,
  subset(ds_T, Condition == "COV")@assays$RNA@data,
  subset(ds_T, Condition == "Sarc")@assays$RNA@data
))

scU_results <- scUnifrac::scUnifrac_multi(
  dataall = dataall,
group = ds_T@meta.data$Condition,
cluster = ds_T@meta.data$seurat_clusters
)

mat <- as.matrix(scU_results$distance)


map <- pheatmap::pheatmap(mat = mat,
                   cluster_rows = T,
                   cluster_cols = T,
                   #breaks = seq(0, max(mat), max(mat)/100)
                   #breaks = c(seq(0, 0.3, 0.00454), seq(0.32, 1, 0.02))
                   breaks = c(seq(0, 0.1, 0.02), seq(0.11, 0.5, 0.0045), seq(0.55, 0.7, 0.17)),
                   color = colorRampPalette(c("tomato", "lightyellow", "steelblue4"))(94)
                   )
ggplot2::ggsave(
  plot = map,
  filename = paste0(plot_path, "BAL_Integrated_T/UniFrac.png"),
  width = 5,
  height = 4
)


######################
#saveRDS(ds_T, file = "/home/alexander/Data/members/Alexander_Leipold/projects/Leo_Rasche/MM_Case_Report/Revision/13_05_22/merged_BALs_T.Rds", compress = FALSE)



Seurat::DimPlot(ds_T, reduction = "mnn_umap", group.by = "seurat_clusters")

Seurat::FeaturePlot(ds_T, reduction = "mnn_umap", features = to_ids["TRGC1"], pt.size = 1, order = T)




























p1 <- plt + ggplot2::geom_point(data = data[data$Condition == "COV", ], ggplot2::aes(x, y, color = density), size = 1.5) + ggplot2::scale_color_gradient2(low = "white", mid = rgb(70, 175, 50, maxColorValue = 255),  high = rgb(15, 45, 15, maxColorValue = 255), midpoint = 0.5)
p2 <- plt + ggplot2::geom_point(data = data[data$Condition == "CP", ], ggplot2::aes(x, y, color = density), size = 1.5) + ggplot2::scale_color_gradient2(low = "white", mid = rgb(70, 175, 50, maxColorValue = 255),  high = rgb(15, 45, 15, maxColorValue = 255), midpoint = 0.5)
p3 <- plt + ggplot2::geom_point(data = data[data$Condition == "HC", ], ggplot2::aes(x, y, color = density), size = 1.5) + ggplot2::scale_color_gradient2(low = "white", mid = rgb(70, 175, 50, maxColorValue = 255),  high = rgb(15, 45, 15, maxColorValue = 255), midpoint = 0.5)
p4 <- plt + ggplot2::geom_point(data = data[data$Condition == "Sarc", ], ggplot2::aes(x, y, color = density), size = 1.5) + ggplot2::scale_color_gradient2(low = "white", mid = rgb(70, 175, 50, maxColorValue = 255),  high = rgb(15, 45, 15, maxColorValue = 255), midpoint = 0.5)

p1 + p2 + p3 + p4















Seurat::DimPlot(ds_T, reduction = "mnn_umap", group.by = "seurat_clusters", pt.size = 2)






Seurat::VlnPlot(ds_T, features = "percent.rib", , group.by = "CP_Ann") +
  ggplot2::geom_violin(draw_quantiles = 0.5, scale = "width", alpha = 0.7)


Seurat::DimPlot(ds_T, reduction = "mnn_umap", group.by = "seurat_clusters", split.by = "batch")

Seurat::FeaturePlot(ds_T, reduction = "mnn_umap", features = to_ids["CD6"], pt.size = 1, order = T)

Seurat::VlnPlot(ds_T, features = to_ids["CD69"])
Seurat::VlnPlot(ds_T, features = "UCell_exTh17", group.by = "CP_Ann")


markers <- Seurat::FindAllMarkers(ds_T, only.pos = T)
markers$gene <- to_genes[markers$gene]
View(markers)

get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

data <- data.frame("cell" = rownames(ds_T@meta.data),
                   "x" = ds_T@reductions$mnn_umap@cell.embeddings[,1],
                   "y" = ds_T@reductions$mnn_umap@cell.embeddings[,2],
                   "Condition" = ds_T@meta.data$Condition)

ggplot2::ggplot(data = data) + ggplot2::geom_point(ggplot2::aes(x, y))

data$density <- NA
data[data$Condition == "COV",]$density <- scales::rescale(get_density(data[data$Condition == "COV",]$x, data[data$Condition == "COV",]$y, n = 100))
data[data$Condition == "CP",]$density <- scales::rescale(get_density(data[data$Condition == "CP",]$x, data[data$Condition == "CP",]$y, n = 100))
data[data$Condition == "HC",]$density <- scales::rescale(get_density(data[data$Condition == "HC",]$x, data[data$Condition == "HC",]$y, n = 100))
data[data$Condition == "Sarc",]$density <- scales::rescale(get_density(data[data$Condition == "Sarc",]$x, data[data$Condition == "Sarc",]$y, n = 100))

plt <- ggplot2::ggplot(data) + ggplot2::geom_point(ggplot2::aes(x, y), color = "gray60", size = 1)
p1 <- plt + ggplot2::geom_point(data = data[data$Condition == "COV", ], ggplot2::aes(x, y, color = density), size = 1.5) + ggplot2::scale_color_gradient2(low = "white", mid = rgb(70, 175, 50, maxColorValue = 255),  high = rgb(15, 45, 15, maxColorValue = 255), midpoint = 0.5)
p2 <- plt + ggplot2::geom_point(data = data[data$Condition == "CP", ], ggplot2::aes(x, y, color = density), size = 1.5) + ggplot2::scale_color_gradient2(low = "white", mid = rgb(70, 175, 50, maxColorValue = 255),  high = rgb(15, 45, 15, maxColorValue = 255), midpoint = 0.5)
p3 <- plt + ggplot2::geom_point(data = data[data$Condition == "HC", ], ggplot2::aes(x, y, color = density), size = 1.5) + ggplot2::scale_color_gradient2(low = "white", mid = rgb(70, 175, 50, maxColorValue = 255),  high = rgb(15, 45, 15, maxColorValue = 255), midpoint = 0.5)
p4 <- plt + ggplot2::geom_point(data = data[data$Condition == "Sarc", ], ggplot2::aes(x, y, color = density), size = 1.5) + ggplot2::scale_color_gradient2(low = "white", mid = rgb(70, 175, 50, maxColorValue = 255),  high = rgb(15, 45, 15, maxColorValue = 255), midpoint = 0.5)

p1 + p2 + p3 + p4

Seurat::FeaturePlot(ds_T, reduction = "mnn_umap", features = to_ids["CD69"], order = T)


genesets <- list(
  exTh17_1 = to_ids[c("CCR6", "DPP4", "CTSH", "LGALS3", "BHLHE40", "PDE4D")],
  exTh17_2 = to_ids[c("TBX21", "IFNG", "CXCR3", "KLRB1", "CCL20", "IL23R",
                      "CCR6", "RORC", "IL4I1", "CSF2", "LGALS3")]
)

scores <- UCell::ScoreSignatures_UCell(ds_T@assays$RNA@data, features = genesets, maxRank = 1500)
head(scores)
ds_T@meta.data$UCell_exTh17_Ray <- scores[,1]
ds_T@meta.data$UCell_exTh17_own <- scores[,2]



CT_ids <- to_ids[unique(c(genes_w_Th17$gene,  as.character(unlist(exth17up))))]

CT_ids <- to_ids[c("CCR6", "DPP4", "CTSH", "LGALS3", "BHLHE40", "PDE4D")]

CT <- list()
CT$Score <- CT_ids


ds_T <- Seurat::AddModuleScore(
  object = ds_T, features = CT, name = names(CT), nbin = 24
)



Seurat::FeaturePlot(ds_T, reduction = "mnn_umap", features = "UCell_exTh17", split.by = "Condition", cols = c("yellow", "darkblue"))

Seurat::VlnPlot(ds_T, idents = c("0"), features = "UCell_exTh17", , group.by = "Condition") +
  ggplot2::geom_violin(draw_quantiles = 0.5, scale = "width", alpha = 0.7)

Seurat::VlnPlot(ds_T, features = "UCell_exTh17_Ray", , group.by = "seurat_clusters") +
  ggplot2::geom_violin(draw_quantiles = 0.5, scale = "width", alpha = 0.7)

Seurat::VlnPlot(ds_T, features = "Score1", , group.by = "seurat_clusters") +
  ggplot2::geom_violin(draw_quantiles = 0.5, scale = "width", alpha = 0.7)

Seurat::VlnPlot(ds_T, features = "Score1", , group.by = "Condition") +
  ggplot2::geom_violin(draw_quantiles = 0.5, scale = "width", alpha = 0.7)

data <- data.frame("Cell" = rownames(ds_T@meta.data),
                   "RORC" = ds_T@assays$RNA@data[to_ids["RORC"],],
                   "TBX21" = ds_T@assays$RNA@data[to_ids["TBX21"],],
                   "CCL20" = ds_T@assays$RNA@data[to_ids["CCL20"],],
                   "IFNG" = ds_T@assays$RNA@data[to_ids["IFNG"],],
                   "CCR6" = ds_T@assays$RNA@data[to_ids["CCR6"],],
                   "CXCR3" = ds_T@assays$RNA@data[to_ids["CXCR3"],])
# get cells that are coexpressing one of the gene-pairs
# RORC/TBX21 -> transcription factors - TF
# CCL20/IFNG -> cytokines - Cyt
# CCR6/CXCR3 -> Receptors - rec

double_pos_TF_cells <- as.character(data$Cell[data$RORC > 0 & data$TBX21 > 0])
double_pos_Cytokine_cells <- as.character(data$Cell[data$CCL20 > 0 & data$IFNG > 0])
double_pos_Receptor_cells <- as.character(data$Cell[data$CCR6 > 0 & data$CXCR3 > 0])

# genuine by coexpression; exact
# add meta.data colums for coexpression of gene pairs - boolean
ds_T@meta.data$genuineRec <- "False"
ds_T@meta.data$genuineRec[rownames(ds_T@meta.data) %in% double_pos_Receptor_cells] <- "True"
ds_T@meta.data$genuineCyt <- "False"
ds_T@meta.data$genuineCyt[rownames(ds_T@meta.data) %in% double_pos_Cytokine_cells] <- "True"
ds_T@meta.data$genuineTF <- "False"
ds_T@meta.data$genuineTF[rownames(ds_T@meta.data) %in% double_pos_TF_cells] <- "True"

ds_T@meta.data$Coexpr <- "No Co-expression"
ds_T@meta.data$Coexpr[ds_T@meta.data$genuineTF == "True"] <- "Single TF"
ds_T@meta.data$Coexpr[ds_T@meta.data$genuineRec == "True"] <- "Single Rec"
ds_T@meta.data$Coexpr[ds_T@meta.data$genuineCyt == "True"] <- "Single Cyt"
ds_T@meta.data$Coexpr[ds_T@meta.data$genuineCyt == "True" &
                              ds_T@meta.data$genuineTF == "True"] <- "Double-TF+Cyt"
ds_T@meta.data$Coexpr[ds_T@meta.data$genuineRec == "True" &
                              ds_T@meta.data$genuineTF == "True"] <- "Double-TF+Rec"
ds_T@meta.data$Coexpr[ds_T@meta.data$genuineRec == "True" &
                              ds_T@meta.data$genuineCyt == "True"] <- "Double-Cyt+Rec"
ds_T@meta.data$Coexpr[ds_T@meta.data$genuineRec == "True" &
                              ds_T@meta.data$genuineCyt == "True" &
                              ds_T@meta.data$genuineTF == "True"] <- "Triple"
ds_T@meta.data$Coexpr <- as.factor(ds_T@meta.data$Coexpr)
ds_T@meta.data$Coexpr <- factor(ds_T@meta.data$Coexpr, levels(ds_T@meta.data$Coexpr)[c(4, 6, 5, 7, 1, 3, 2, 8)])




td <- table(ds_T@meta.data$Condition[ds_T@meta.data$seurat_clusters == "1"], ds_T@meta.data$Coexpr[ds_T@meta.data$seurat_clusters == "1"])
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
colnames(data) <- c("Cluster", "CoexprType", "proportion")

plot <- plot <- ggplot2::ggplot(
  data    = data,
  mapping = ggplot2::aes(
    x     = Cluster,
    y     = proportion,
    col   = CoexprType,
    fill  = CoexprType
  )
) +
  ggplot2::geom_col(position = "stack") +
  ggplot2::theme(
    panel.background = ggplot2::element_blank(),
    panel.border     = ggplot2::element_rect(size = 2, fill = NA),
    axis.text        = ggplot2::element_text(size = 15),
    axis.text.x      = ggplot2::element_text(), 
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



###################################################################################################

ds_T_CP_HC <- subset(ds_T, subset = Condition %in% c("CP", "HC", "Sarc"))


ds_T_CP_HC <- Seurat::NormalizeData(
  object               = ds_T_CP_HC,
  assay                = "RNA",
  normalization.method = "LogNormalize",
  scale.factor         = 10000
)
ds_T_CP_HC <- Seurat::FindVariableFeatures(
  object           = ds_T_CP_HC,
  selection.method = "vst", 
  nfeatures        = 3000
)


set.seed(seed = seed)
ds_T_CP_HC@reductions[["mnn"]] <- Seurat::CreateDimReducObject(
  embeddings = batchelor::fastMNN(
    ds_T_CP_HC@assays$RNA@data[ds_T_CP_HC@assays$RNA@var.features, ],
    batch = ds_T_CP_HC@meta.data$Condition,
    d     = 55
  )@int_colData$reducedDims$corrected,
  key        = "MNN_",
  assay      = "RNA"
)


# UMAP with MNN-corrected PCA
ds_T_CP_HC <- Seurat::RunUMAP(
  object         = ds_T_CP_HC,
  dims           = 1:55,
  seed.use       = seed,
  reduction      = "mnn",
  reduction.name = "mnn_umap")


# Clustering with Louvain algorithm based on SNN-graph
ds_T_CP_HC <- Seurat::FindNeighbors(
  object    = ds_T_CP_HC,
  dims      = 1:55,
  reduction = "mnn"
)

ds_T_CP_HC <- Seurat::FindClusters(
  object     = ds_T_CP_HC,
  resolution = 0.2
)


Seurat::DimPlot(ds_T_CP_HC, group.by = "seurat_clusters", reduction = "mnn_umap", split.by = "Condition")


ds_T_CP_HC <- Seurat::SetIdent(ds_T_CP_HC, value = ds_T_CP_HC@meta.data$Condition)

mrkrs <- Seurat::FindMarkers(subset(ds_T_CP_HC, subset = seurat_clusters == "0"), ident.1 = "CP", ident.2 = "HC", only.pos = TRUE)
mrkrs$gene <- to_genes[rownames(mrkrs)]
View(mrkrs)

Seurat::FeaturePlot(ds_T_CP_HC, reduction = "mnn_umap", features = to_ids["IL26"])

genesets <- list(
  exTh17 = to_ids[c("CCR6", "DPP4", "CTSH", "LGALS3", "BHLHE40", "PDE4D")]
)

scores <- UCell::ScoreSignatures_UCell(ds_T_CP_HC@assays$RNA@data, features=genesets)
head(scores)
ds_T_CP_HC@meta.data$UCell_exTh17 <- scores[,1]

Seurat::VlnPlot(ds_T_CP_HC, features = "UCell_exTh17", group.by = "Condition", split.by = "seurat_clusters")




data <- data.frame("Cell" = rownames(ds_T_CP_HC@meta.data),
                   "RORC" = ds_T_CP_HC@assays$RNA@data[to_ids["RORC"],],
                   "TBX21" = ds_T_CP_HC@assays$RNA@data[to_ids["TBX21"],],
                   "CCL20" = ds_T_CP_HC@assays$RNA@data[to_ids["CCL20"],],
                   "IFNG" = ds_T_CP_HC@assays$RNA@data[to_ids["IFNG"],],
                   "CCR6" = ds_T_CP_HC@assays$RNA@data[to_ids["CCR6"],],
                   "CXCR3" = ds_T_CP_HC@assays$RNA@data[to_ids["CXCR3"],])
# get cells that are coexpressing one of the gene-pairs
# RORC/TBX21 -> transcription factors - TF
# CCL20/IFNG -> cytokines - Cyt
# CCR6/CXCR3 -> Receptors - rec

double_pos_TF_cells <- as.character(data$Cell[data$RORC > 0 & data$TBX21 > 0])
double_pos_Cytokine_cells <- as.character(data$Cell[data$CCL20 > 0 & data$IFNG > 0])
double_pos_Receptor_cells <- as.character(data$Cell[data$CCR6 > 0 & data$CXCR3 > 0])

# genuine by coexpression; exact
# add meta.data colums for coexpression of gene pairs - boolean
ds_T_CP_HC@meta.data$genuineRec <- "False"
ds_T_CP_HC@meta.data$genuineRec[rownames(ds_T_CP_HC@meta.data) %in% double_pos_Receptor_cells] <- "True"
ds_T_CP_HC@meta.data$genuineCyt <- "False"
ds_T_CP_HC@meta.data$genuineCyt[rownames(ds_T_CP_HC@meta.data) %in% double_pos_Cytokine_cells] <- "True"
ds_T_CP_HC@meta.data$genuineTF <- "False"
ds_T_CP_HC@meta.data$genuineTF[rownames(ds_T_CP_HC@meta.data) %in% double_pos_TF_cells] <- "True"

ds_T_CP_HC@meta.data$Coexpr <- "No Co-expression"
ds_T_CP_HC@meta.data$Coexpr[ds_T_CP_HC@meta.data$genuineTF == "True"] <- "Single TF"
ds_T_CP_HC@meta.data$Coexpr[ds_T_CP_HC@meta.data$genuineRec == "True"] <- "Single Rec"
ds_T_CP_HC@meta.data$Coexpr[ds_T_CP_HC@meta.data$genuineCyt == "True"] <- "Single Cyt"
ds_T_CP_HC@meta.data$Coexpr[ds_T_CP_HC@meta.data$genuineCyt == "True" &
                      ds_T_CP_HC@meta.data$genuineTF == "True"] <- "Double-TF+Cyt"
ds_T_CP_HC@meta.data$Coexpr[ds_T_CP_HC@meta.data$genuineRec == "True" &
                      ds_T_CP_HC@meta.data$genuineTF == "True"] <- "Double-TF+Rec"
ds_T_CP_HC@meta.data$Coexpr[ds_T_CP_HC@meta.data$genuineRec == "True" &
                      ds_T_CP_HC@meta.data$genuineCyt == "True"] <- "Double-Cyt+Rec"
ds_T_CP_HC@meta.data$Coexpr[ds_T_CP_HC@meta.data$genuineRec == "True" &
                      ds_T_CP_HC@meta.data$genuineCyt == "True" &
                      ds_T_CP_HC@meta.data$genuineTF == "True"] <- "Triple"
ds_T_CP_HC@meta.data$Coexpr <- as.factor(ds_T_CP_HC@meta.data$Coexpr)
ds_T_CP_HC@meta.data$Coexpr <- factor(ds_T_CP_HC@meta.data$Coexpr, levels(ds_T_CP_HC@meta.data$Coexpr)[c(4, 6, 5, 7, 1, 3, 2, 8)])


td <- table(ds_T_CP_HC@meta.data$Condition[ds_T_CP_HC@meta.data$seurat_clusters == "0"], ds_T_CP_HC@meta.data$Coexpr[ds_T_CP_HC@meta.data$seurat_clusters == "0"])
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
colnames(data) <- c("Cluster", "CoexprType", "proportion")

plot <- plot <- ggplot2::ggplot(
  data    = data,
  mapping = ggplot2::aes(
    x     = Cluster,
    y     = proportion,
    col   = CoexprType,
    fill  = CoexprType
  )
) +
  ggplot2::geom_col(position = "stack") +
  ggplot2::theme(
    panel.background = ggplot2::element_blank(),
    panel.border     = ggplot2::element_rect(size = 2, fill = NA),
    axis.text        = ggplot2::element_text(size = 15),
    axis.text.x      = ggplot2::element_text(), 
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



###################################

ds.list <- Seurat::SplitObject(ds_T, split.by = "Origin")
ds.list <- lapply(X = ds.list, FUN = function(x) {
  x <- Seurat::NormalizeData(x)
  x <- Seurat::FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)
})

ds_T@assays$RNA@var.features <- Seurat::SelectIntegrationFeatures(object.list = ds.list, nfeatures = 2000)




CT_CPBAL <- ds@meta.data$seurat_clusters
names(CT_CPBAL) <- rownames(ds@meta.data)
CT_CPBAL <- CT_CPBAL[colnames(ds[,colnames(ds) %in% colnames(subset(ds_T, subset = Condition == "CP"))])]

ds_T@meta.data$CP_Ann <- "other"
ds_T@meta.data[names(CT_CPBAL),]$CP_Ann <- as.character(as.character(CT_CPBAL))

Seurat::DimPlot(subset(ds_T, subset = Condition == "CP"), reduction = "mnn_umap", group.by = "CP_Ann", pt.size = 2)





genexl <- readxl::read_excel("/home/alexander/Desktop/1-s2.0-S1074761322001753-mmc5.xlsx", sheet = 14)
gene_ids_pos <- to_ids[genexl$Gene_ID[genexl$avg_logFC > 0 & genexl$p_val_adj < 1e-100 & genexl$avg_logFC > 0.5]]
gene_ids_neg <- to_ids[genexl$Gene_ID[genexl$avg_logFC < 0 & genexl$p_val_adj < 1e-100 & genexl$avg_logFC < -0.5]]
gene_ids_pos <- gene_ids_pos[!is.na(gene_ids_pos)]
gene_ids_neg <- gene_ids_neg[!is.na(gene_ids_neg)]
gene_ids_pos <- gene_ids_pos[gene_ids_pos %in% ds@assays$RNA@var.features]
gene_ids_neg <- gene_ids_neg[gene_ids_neg %in% ds@assays$RNA@var.features]

gene_ids <- c(paste0(gene_ids_pos, "+"), paste0(gene_ids_neg, "-"))



gene_ids <- to_ids[c("RORA", "RORC", "RBPJ", "BHLHE40", "FURIN", "COTL1", "CCL20", "CCR6",
                     "IL23R", "CXCR3", "IFNG", "TNF")]

CT <- list()
CT$exTh17ScoreMac <- gene_ids

ds_T <- Seurat::AddModuleScore(
  object = ds_T, features = CT, name = names(CT), nbin = 24
)

Seurat::VlnPlot(ds_T, features = "exTh17Score1", group.by = "seurat_clusters") +
  ggplot2::geom_violin(draw_quantiles = 0.5, scale = "width", alpha = 0.75)

Seurat::DimPlot(ds_T, reduction = "mnn_umap", group.by = "ds_T$RNA_snn_res.0.4")






Seurat::VlnPlot(ds_T, features = "nCount_RNA", group.by = "Condition") +
  ggplot2::geom_violin(draw_quantiles = 0.5, scale= "width", alpha = 0.7)











############################
clust.tree <- function(object, method = "pca", group = "seurat_clusters",
                       dim_red_use = "mnn", pcs_use = NULL, genes_use = NULL) {
  require(Seurat)
  require(ape)
  
  if(method == "pca") {
    data_use <- data.frame(row.names = c(1:pcs_use))
    
    for (i in levels(object@meta.data[[group]])) {
      pca_tmp <- object@reductions[[dim_red_use]]@cell.embeddings[ds_T@meta.data$seurat_clusters == i, 1:pcs_use]
      data_tmp <- apply(pca_tmp, 2, mean)
      data_use <- cbind(data_use, data_tmp)
    }
    colnames(data_use) <- paste0("c", levels(object@meta.data[[group]]))
    
    weights = object@reductions$pca@stdev[1:pcs_use]
    data_dist <- Seurat::CustomDistance(my.mat = data_use, my.function = weighted.euclidean.distance, w = weights)
    
    data_tree <- hclust(data_dist)
    data_tree_phylo <- ape::as.phylo(data_tree)
    return(list(hclust = data_tree, plot = plot.phylo(data_tree_phylo, direction = "downwards")))
  }
}

weighted.euclidean.distance <- function(x, y, w) {
  v.dist <- sum(sqrt(x = w * (x - y) ^ 2))
  return(v.dist)
}

tree <- clust.tree(ds_T, method = "pca", group = "seurat_clusters", pcs_use = 55)



sc.unifrac.multi <- function(object, group_by, clusts_name, tree, perm_iters, n_cores){
  require(doRNG)
  require(GUniFrac)
  
  groups <- setNames(factor(object@meta.data[[group_by]]), colnames(object))
  group_names = levels(groups)
  
  if(!is.factor(groups)) stop("\"group_by\" must be a factor")
  if(length(levels(groups)) == 2) stop("only two levels in the groups - use sc.unifrac.pair")
  
  clusts <- setNames(object@meta.data[[clusts_name]], colnames(object))
  num_clusts <- max(as.numeric(levels(clusts)))
  
  count_table <- table(groups, clusts)
  colnames(count_table) <- paste0("c", levels(clusts))
  
  tree_phylo <- as.phylo(tree$hclust)
  
  unifracs <- GUniFrac(count_table, tree_phylo, alpha = c(0, 0.5, 1))$unifracs
  dist_obs <- unifracs[, , "d_1"]
  rownames(dist_obs) <- colnames(dist_obs) <- rownames(count_table)
 # return(dist_obs)
  
  cl <- parallel::makeCluster(n_cores, outfile = "")
  doParallel::registerDoParallel(cl)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  
  perm_scores <- foreach(i = 1:perm_iters, .packages = c('scUnifrac', "GUniFrac")) %dorng% {
    groups_perm <- factor(setNames(sample(as.character(groups)), names(groups)),
                          levels = levels(groups))
    
    count_table_perm <- table(groups_perm, clusts)
    colnames(count_table_perm) <- paste0("c", levels(clusts))
    
    unifracs_perm <- GUniFrac(count_table_perm, tree_phylo, alpha = c(0, 0.5, 1))$unifracs
    
    return(unifracs_perm[, , "d_1"])
  }
  
  perm_scores_array <- array(unlist(perm_scores), dim = c(length(levels(groups)), length(levels(groups)), perm_iters))
  p_vals <- matrix(apply(apply(perm_scores_array, 3, ">=", dist_obs), 1, sum), nrow = length(levels(groups)), byrow = F)/perm_iters
  rownames(p_vals) <- colnames(p_vals) <- rownames(count_table)
  
  return(list(distance = dist_obs, count_table = count_table, pvalue = p_vals))
}

reslts <- sc.unifrac.multi(object = ds_T, group_by = "Condition", clusts_name = "seurat_clusters", tree = tree, perm_iters = 2, n_cores = 1)
reslts$distance


reslts$distance <- reslts$distance[c(2, 4, 3, 1), c(2, 4, 3, 1)]
mat <- as.matrix(reslts$distance)


map <- pheatmap::pheatmap(mat = mat,
                          cluster_rows = F,
                          cluster_cols = F,
                          #breaks = seq(0, max(mat), max(mat)/100)
                          #breaks = c(seq(0, 0.3, 0.00454), seq(0.32, 1, 0.02)),
                          #breaks = c(seq(0, 0.1, 0.02), seq(0.11, 0.5, 0.0045), seq(0.55, 0.7, 0.17)),
                          #breaks = c(seq(0, 0.04, 0.01), seq(0.05, 0.35, 0.001), seq(0.352, 0.702, 0.01)),
                          color = colorRampPalette(c("tomato", "lightyellow", "steelblue4"))(94)#342)
)
ggplot2::ggsave(
  plot = map,
  filename = paste0(plot_path, "BAL_Integrated_T/UniFrac.png"),
  width = 4.7,
  height = 4
)

index <- c("CART_BAL", "COVID-19", "Sarc", "HC")
names(index) <- c("CP", "COV", "Sarc", "HC")
ds_T@meta.data$Cond <- index[ds_T@meta.data$Condition]
