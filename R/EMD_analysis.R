# Script for analysis of EMD sample
# Alexander Leipold

# set seed for reproducibility
seed <- 1992

# set path for plots
plot_path <- "/path/for/plots/"

# calculate percentage of mitochondrial mRNA Counts per cell
mt_counts <- Matrix::colSums(ds@assays$RNA@counts[
  ds@misc$features$ENSEMBL[which(stringr::str_detect(features$SYMBOL, "^MT-"))]
  , ]) 
ds@meta.data$percent.mt <- round(
  mt_counts / ds@meta.data$nCount_RNA * 100, 1
)


# Quality Thresholds
percent.mt_min <- 0.2
percent.mt_max <- 12.5
nFeature_RNA_min <- 500
nFeature_RNA_max <- 5000
nCount_RNA_min <- 1000
nCount_RNA_max <- 50000


# Set cell (column) filter
# Based on thersholds for QC metrics
cells <- rownames(ds@meta.data)[which(
  ds@meta.data$nFeature_RNA > nFeature_RNA_min &
    ds@meta.data$nFeature_RNA < nFeature_RNA_max &
    ds@meta.data$nCount_RNA > nCount_RNA_min &
    ds@meta.data$nCount_RNA < nCount_RNA_max &
    ds@meta.data$percent.mt > percent.mt_min & 
    ds@meta.data$percent.mt < percent.mt_max
)]

# apply filtering
ds <- subset(ds, cells = cells)

# basic Seurat workflow
# Normalize mRNA Counts and Find highly variable genes
ds <- Seurat::NormalizeData(
  object               = ds,
  normalization.method = "LogNormalize",
  scale.factor         = 10000
)

ds <- Seurat::FindVariableFeatures(
  object           = ds,
  selection.method = "vst", 
  nfeatures        = 1000
)                   

# Scale data 
ds <- Seurat::ScaleData(ds)

# dimensionality reduction using PCA and UMAP
ds <- Seurat::RunPCA(ds, npcs = 50)

ds <- Seurat::RunUMAP(
  object    = ds,
  dims      = 1:10,
  reduction = "pca",
  seed.use  = seed
)

# Clustering with Louvain algorithm based on SNN-graph
ds <- Seurat::FindNeighbors(
  object    = ds,
  dims      = 1:10,
  reduction = "pca"
)


# Plot Celltype
CT.colors <- c(
  "MM" = "indianred3",
  "T/NK" = "skyblue2",
  "Mono/Macro" = "mediumorchid2",
  "Fibroblast" = "springgreen3"
)

# Plot Subsets
embedding <- "umap"
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
    size = 0.9
  ) +
  ggplot2::scale_color_manual(values = CT.colors) +
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
  filename = paste0(plot_path, "CellType.png"),
  width = 6,
  height = 5
)

# Plot expression of interesting genes
genes <- c(
  "CXCR4", "TNFRSF17", "MKI67", "CCND2",
  "EPCAM", "SFN", "KRT8", "KRT18",
  ) 
pltlist = list()
for (i in genes) {
  embedding <- "umap"
  data <- tidyr::as_tibble(
    x = ds@reductions[[embedding]]@cell.embeddings
  )
  colnames(data) <- c("x", "y")
  category <- i
  data[["value"]] <- ds@assays$RNA@data[to_ids[category],]
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
      size = 0.6
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


# UMAPxgene plots Figure2
voidplt <- ggplot2::ggplot() + ggplot2::theme_void()

pltlistFig <- list(plot_CXCR4, plot_TNFRSF17, plot_MKI67, plot_CCND2,
                    voidplt, voidplt, voidplt, voidplt,
                    plot_EPCAM, plot_SFN, plot_KRT8, plot_KRT18)
plot <- ggpubr::ggarrange(plotlist = pltlistFig,
                          nrow = 3, ncol = 4,
                          heights = c(1, -0.45, 1))
ggplot2::ggsave(
  plot = plot,
  filename = paste0(plot_path, "Genes.png"),
  width = 10,
  height = 10
)

########################################################################

# Script for comparison of EMD and IMD MM transcriptomes
# Alexander Leipold

# ds_IMD is a Seurat object containing quality filtered (see Supplementary table 1) IMD samples
ds <- merge(x = ds, y = ds_IMD)

seed = 1992

# basic Seurat workflow
# Normalize mRNA Counts and Find highly variable genes
ds <- Seurat::NormalizeData(
  object               = ds,
  normalization.method = "LogNormalize",
  scale.factor         = 10000
)

ds <- Seurat::FindVariableFeatures(
  object           = ds,
  selection.method = "vst", 
  nfeatures        = 1000
)                   

# Scale data 
ds <- Seurat::ScaleData(ds)

# dimensionality reduction using PCA and UMAP
ds <- Seurat::RunPCA(ds, npcs = 50)

ds <- Seurat::RunUMAP(
  object    = ds,
  dims      = 1:30,
  reduction = "pca",
  seed.use  = seed
)

ds <- Seurat::SetIdent(ds, value = ds@meta.data$MMType)
markers <- Seurat::FindAllMarkers(ds, only.pos = T)
markers$gene <- to_genes[markers$gene]
View(markers)
markers <- markers[,c(6, 7, 2, 3, 4, 1, 5)]
markers$ENSEMBL_ID <- to_ids[markers$gene]
# Save as csv
write.csv(markers, file = "EMDIMD_DE.csv", row.names = FALSE)

# Heatmap of selected genes
# genes were selected due to funtions that are known to play roles in 
# extramedullary niche homing
data <- data.frame(
  "Sample" = d@meta.data$Sample,
  "EPCAM" = ds@assays$RNA@data[to_ids["EPCAM"],],
  "SFN"   = ds@assays$RNA@data[to_ids["SFN"],],
  "KRT8"   = ds@assays$RNA@data[to_ids["KRT8"],],
  "KRT18"   = ds@assays$RNA@data[to_ids["KRT18"],],
  "ADGRE2"   = ds@assays$RNA@data[to_ids["ADGRE2"],],
  "ITGB7"   = ds@assays$RNA@data[to_ids["ITGB7"],],
  "VEGFA"   = ds@assays$RNA@data[to_ids["VEGFA"],],
  "PTN"   = ds@assays$RNA@data[to_ids["PTN"],],
  "ANG"   = ds@assays$RNA@data[to_ids["ANG"],],
  "TGFA"   = ds@assays$RNA@data[to_ids["TGFA"],]
  )

library(dplyr)
data <- data %>% group_by(Sample) %>% summarise(across(everything(), list(mean)))
colnames(data) <- c("Sample", gsub("_1", "", colnames(data[, 2:length(colnames(data))])))
data <- as.matrix(data[, 2:length(colnames(data))])
data <- scale(data)
data[data > 2] <- 2
data[data < -2] <- -2

rann <- data.frame("Cluster" = levels(ds_T@meta.data$T_subset))
rownames(data) <- levels(ds_EMDIMD@meta.data$Sample)
rownames(rann) <- rownames(data)
data <- t(data)
map <- pheatmap::pheatmap(mat = data,
                          color = colorRampPalette(rev(RColorBrewer::brewer.pal(n=7, name = "RdBu")))(500)[1:500],
                          scale = "none",
                          cluster_rows = F,
                          cluster_cols = F,
                          annotation_names_row = TRUE
)

filename <- paste0(
  plot_path,
  "Heatmap.png",
  ".png"
)

ggplot2::ggsave(
  filename = filename,
  plot     = map,
  width    = 10,
  height   = 5
)


