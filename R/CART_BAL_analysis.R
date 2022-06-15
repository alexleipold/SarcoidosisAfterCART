# Script for the analysis of BAL fluid scRNA-seq data of case patient - CP
# technical duplicate

# clean environment
rm(list = ls())
gc()

# set seed for reproducibility
seed <- 1999

# set path for plots
plot_path <- "/home/alexander/Data/members/Alexander_Leipold/projects/Leo_Rasche/MM_Case_Report/Revision/13_05_22/plots/"

datasets <- c("CP_BAL/A8_3", "CP_BAL/A9_3")

# Load datasets_samples and add meta.data that we need later

for (dataset in datasets) {
  from <- paste0("/home/alexander/Data/members/Alexander_Leipold/projects/Leo_Rasche/MM_Case_Report/Revision/13_05_22/datasets/", dataset)

  features <- readr::read_tsv(
    file = gzfile(paste0(from, "/filtered_feature_bc_matrix/features.tsv.gz")),
    col_names = c("ENSEMBL", "SYMBOL", "TYPE")
  )
  barcodes <- readr::read_lines(
    file = gzfile(paste0(from, "/filtered_feature_bc_matrix/barcodes.tsv.gz"))
  )
  ds <- Matrix::readMM(
    file = paste0(from, "/SoupX_matrix.mtx")
  )
  
  rownames(ds) <- features$ENSEMBL
  colnames(ds) <- barcodes
  
  # Create Seurat object from matrix
  ds <- Seurat::CreateSeuratObject(
    counts = ds
  )
  # and store gene/ensID conversion table
  ds@misc$features <- features
  ds@meta.data$Origin <- dataset
  
  # assign name related to dataset
  assign(paste0("ds_", dataset), ds)
  
}

# merge datasets
ds <- merge(x = `ds_CP_BAL/A8_3`, y = `ds_CP_BAL/A9_3`)

# add meta.data technical replicate data
index <- c("CP_BAL_replicate1", "CP_BAL_replicate2")
names(index) <- c("CP_BAL/A8_3", "CP_BAL/A9_3")
ds@meta.data$technical_replicate <- index[ds@meta.data$Origin]

# add Ensembl/symbol reference to ds@misc
ds@misc$features <- `ds_CP_BAL/A8_3`@misc$features

# create gene symbol-EnsemblID lookup table
to_ids <- ds@misc$features$ENSEMBL
names(to_ids) <- ds@misc$features$SYMBOL

to_genes <- ds@misc$features$SYMBOL
names(to_genes) <- ds@misc$features$ENSEMBL

# calculate percentage of mitochondrial mRNA Counts per cell
mt_counts <- Matrix::colSums(ds@assays$RNA@counts[
  ds@misc$features$ENSEMBL[which(stringr::str_detect(features$SYMBOL, "^MT-"))]
  , ]) 
ds@meta.data$percent.mt <- round(
  mt_counts / ds@meta.data$nCount_RNA * 100, 1
)

# Quality Thresholds
percent.mt_min <- 0.2
percent.mt_max <- 10
nFeature_RNA_min <- 750
nFeature_RNA_max <- 7500
nCount_RNA_min <- 1300
nCount_RNA_max <- 90000

# Data Quality plots
# UMI Counts
quality_metric <- "nCount_RNA"
data <- tidyr::as_tibble(t(as.matrix(
  ds@meta.data[[quality_metric]]
)))
data <- data.frame("Value" = ds@meta.data[[quality_metric]])
data$replicate <- ds@meta.data$technical_replicate
plot <- ggplot2::ggplot(
  data    = data,
  mapping = ggplot2::aes(
    x     = replicate,
    y     = Value,
    fill  = replicate
  )
) +
  ggplot2::geom_point(
    position = "jitter", size = 0.2, color = "black", 
  ) +
  ggplot2::geom_violin(
    scale = "width", size = 0.5, alpha = 0.75, draw_quantiles = 0.5,
    trim = TRUE
  ) +
  ggplot2::scale_color_distiller(palette = "Reds", direction = 1) +
  ggplot2::theme_classic(base_size = 15) +
  ggplot2::scale_y_log10(breaks = c(1000, 10000, 100000), labels = c("1k", "10k", "100k")) +
  ggplot2::geom_hline(yintercept = nCount_RNA_max, colour = "red", size = 1.25) +
  ggplot2::geom_hline(yintercept = nCount_RNA_min, colour = "red", size = 1.25) +
  ggplot2::theme(
    legend.position  = "",
    axis.title       = ggplot2::element_blank(),
    axis.text.x      = ggplot2::element_text(angle = 45, hjust = 1),
    strip.background = ggplot2::element_blank()
  )
ggplot2::ggsave(
  filename = paste0(plot_path, "DataQuality/QC_Counts.png"),
  width = 4,
  height = 4.4
)

# Genes per cell
quality_metric <- "nFeature_RNA"
data <- tidyr::as_tibble(t(as.matrix(
  ds@meta.data[[quality_metric]]
)))
data <- data.frame("Value" = ds@meta.data[[quality_metric]])
data$replicate <- ds@meta.data$technical_replicate
plot <- ggplot2::ggplot(
  data    = data,
  mapping = ggplot2::aes(
    x     = replicate,
    y     = Value,
    fill  = replicate
  )
) +
  ggplot2::geom_point(
    position = "jitter", size = 0.3, color = "black", 
  ) +
  ggplot2::geom_violin(
    scale = "width", size = 0.5, alpha = 0.75, draw_quantiles = 0.5,
    trim = TRUE
  ) +
  ggplot2::scale_color_distiller(palette = "Reds", direction = 1) +
  ggplot2::theme_classic(base_size = 15) +
  ggplot2::scale_y_log10(breaks = c(100, 1000, 10000), labels = c("100", "1,000", "10,000")) +
  ggplot2::geom_hline(yintercept = nFeature_RNA_max, colour = "red", size = 1.25) +
  ggplot2::geom_hline(yintercept = nFeature_RNA_min, colour = "red", size = 1.25) +
  ggplot2::theme(
    legend.position  = "",
    axis.title       = ggplot2::element_blank(),
    axis.text.x      = ggplot2::element_text(angle = 45, hjust = 1),
    strip.background = ggplot2::element_blank()
  )
ggplot2::ggsave(
  filename = paste0(plot_path, "DataQuality/QC_Features.png"),
  width = 4,
  height = 4.4
)

# mitochondrial count fraction
quality_metric <- "percent.mt"
data <- tidyr::as_tibble(t(as.matrix(
  ds@meta.data[[quality_metric]]
)))
data <- data.frame("Value" = ds@meta.data[[quality_metric]])
data$replicate <- ds@meta.data$technical_replicate
plot <- ggplot2::ggplot(
  data    = data,
  mapping = ggplot2::aes(
    x     = replicate,
    y     = Value,
    fill  = replicate
  )
) +
  ggplot2::geom_point(
    position = "jitter", size = 0.3, color = "black", 
  ) +
  ggplot2::geom_violin(
    scale = "width", size = 0.5, alpha = 0.75, draw_quantiles = 0.5,
    trim = TRUE
  ) +
  ggplot2::scale_color_distiller(palette = "Reds", direction = 1) +
  ggplot2::theme_classic(base_size = 15) +
  ggplot2::geom_hline(yintercept = percent.mt_max, colour = "red", size = 1.25) +
  ggplot2::geom_hline(yintercept = percent.mt_min, colour = "red", size = 1.25) +
  ggplot2::theme(
    legend.position  = "",
    axis.title       = ggplot2::element_blank(),
    axis.text.x      = ggplot2::element_text(angle = 45, hjust = 1),
    strip.background = ggplot2::element_blank()
  )
ggplot2::ggsave(
  filename = paste0(plot_path, "DataQuality/QC_pctmito.png"),
  width = 4,
  height = 4.4
)

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

# identify epithelial cells based on EPCAM expression
cells <- colnames(ds[, ds@assays$RNA@data[to_ids["EPCAM"],] == 0])
# apply filtering of epithelial cells
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
  nfeatures        = 5000
)                   

# Scale data 
ds <- Seurat::ScaleData(ds)

# dimensionality reduction using PCA and UMAP
ds <- Seurat::RunPCA(ds, npcs = 50)

ds <- Seurat::RunUMAP(
  object    = ds,
  dims      = 1:45,
  reduction = "pca",
  seed.use  = seed
)

# Clustering with Louvain algorithm based on SNN-graph
ds <- Seurat::FindNeighbors(
  object    = ds,
  dims      = 1:45,
  reduction = "pca"
)

ds <- Seurat::FindClusters(
  object     = ds,
  resolution = 1.8
)

# Plot Replicates
embedding <- "umap"
data <- tidyr::as_tibble(
  x = ds@reductions[[embedding]]@cell.embeddings
)
colnames(data) <- c("x", "y")
category <- "technical_replicate"
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
  filename = paste0(plot_path, "DataQuality/Replicate_CP_BAL.png"),
  width = 7.5,
  height = 5
)

# Plot Clusters
embedding <- "umap"
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
  filename = paste0(plot_path, "DataQuality/Cluster_CP_BAL.png"),
  width = 6.5,
  height = 5
)

# Differential expression analysis Clusters
markers <- Seurat::FindAllMarkers(ds, only.pos = TRUE)
markers$gene <- to_genes[markers$gene]
markers <- markers[,c(6, 7, 2, 3, 4, 1, 5)]
markers$ENSEMBL_ID <- to_ids[markers$gene]
# Save as csv
write.csv(markers, file = "/home/alexander/Data/members/Alexander_Leipold/projects/Leo_Rasche/MM_Case_Report/Revision/13_05_22/tables/CP_BAL_Cluster_DE.csv", row.names = FALSE)

# Dotplot of marker genes
# marker genes
marker_genes <- c("CD68", "CD163", "CD14", "MARCO", "FABP4", # Mp
                  "FCN1", "VCAN", # Mono/Mp
                  "MKI67", "TOP2A", #Prolif
                  "FCER1A", "CD1C", "CD1E", "FLT3", # mDC
                  "BASP1", "LGALS2", "CLEC9A", #pDC
                  "FCGR3B", "CXCL8", "CSF3R", #Neutrophils
                  "TPSAB1", "MS4A2", "GATA2", #Mast
                  "CD3D", "CD3E", "IL32", # pan T-cell
                  "FOXP3", "IL2RA", # Treg
                  "CD4", "CD40LG", #CD4
                  "CD8A", "CD8B", # CD8
                  "TRDC", "KLRC1", "XCL1" # NK
)
genes <- rev(marker_genes)
data <- t(as.matrix(ds@assays$RNA@data[to_ids[genes], ]))
data <- tidyr::as_tibble(data)
colnames(data) <- genes
data$Cluster <- factor(ds@meta.data$seurat_clusters)
data$Cluster <- factor(data$Cluster, levels(data$Cluster)[c(3,4,9,10,11,12,16, # AMp
                                                            5, # Mono/Mp
                                                            15, 21, # Prolif Mp
                                                            18, # mDC
                                                            20, # pDC
                                                            22, # Neutro
                                                            24, # Mast
                                                            14, # Treg
                                                            1, 2, 6, 13, 17, # CD4 T
                                                            23, # Prolif T
                                                            7, 8, # CD8 T
                                                            19 # NK
)])
data$Cluster <- factor(data$Cluster, levels(data$Cluster)[seq(24,1)])
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
  filename = paste0(plot_path, "DataQuality/DP.png"),
  width = 9,
  height = 6
)

# Annotate Clusters
cluster.annotation <- c(
  "0" = "T CD4",
  "1" = "T CD4",
  "2" = "AMp",
  "3" = "AMp",
  "4" = "Mono/Mp",
  "5" = "T CD4",
  "6" = "T CD8",
  "7" = "T CD8", 
  "8" = "AMp",
  "9" = "AMp",
  "10" = "AMp",
  "11" = "AMp",
  "12" = "T CD4",
  "13" = "Treg",
  "14" = "AMp Prolif",
  "15" = "AMp",
  "16" = "T CD4",
  "17" = "mDC",
  "18" = "NK",
  "19" = "pDC",
  "20" = "AMp Prolif",
  "21" = "Neutro",
  "22" = "T Prolif",
  "23" = "Mast"
)
ds@meta.data$Celltype <- factor(
  cluster.annotation[as.character(ds$seurat_clusters)],
  unique(cluster.annotation)[c(2, 3, 6, 7, 9, 10, 12, 5, 1, 11, 4, 8)]
)

# Differential expression analysis Celltype
ds <- Seurat::SetIdent(ds, value = ds@meta.data$Celltype)
markers <- Seurat::FindAllMarkers(ds, only.pos = TRUE)
markers$gene <- to_genes[markers$gene]
markers <- markers[,c(6, 7, 2, 3, 4, 1, 5)]
markers$ENSEMBL_ID <- to_ids[markers$gene]
# Save as csv
write.csv(markers, file = "/home/alexander/Data/members/Alexander_Leipold/projects/Leo_Rasche/MM_Case_Report/Revision/13_05_22/tables/CP_BAL_Celltypes_DE.csv", row.names = FALSE)
ds <- Seurat::SetIdent(ds, value = ds@meta.data$seurat_clusters)

# Plot Celltype
CT_Colors <- c(
  "AMp"         = "mediumorchid1",
  "Mono/Mp"     = "mediumpurple",
  "AMp Prolif"  = "mediumorchid3",
  "mDC"         = "plum1",
  "pDC"         = "plum3",
  "Neutro"  = "yellow4",
  "Mast"        = "turquoise",
  "Treg"        = "skyblue4",
  "T CD4"       = "skyblue1", 
  "T CD8"       = "steelblue3",
  "NK"          = "slateblue",
  "T Prolif"    = "cyan4"
)
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
      reverse = TRUE
    )) +
  ggplot2::scale_color_manual(
    values = CT_Colors[rev(levels(data$col))]
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
  filename = paste0(plot_path, "DataQuality/CT_CP_BAL.png"),
  width = 6.5,
  height = 5
)

# Plot quality according to Celltype
# UMI Counts per cell
quality_metric <- "nCount_RNA"
data <- tidyr::as_tibble(t(as.matrix(
  ds@meta.data[[quality_metric]]
)))
data <- data.frame("Value" = ds@meta.data[[quality_metric]])
data$col <- ds@meta.data$Celltype
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
  ggplot2::geom_violin(
    scale = "width", size = 0.5, alpha = 0.75, draw_quantiles = 0.5,
    trim = TRUE
  ) +
  ggplot2::scale_fill_manual(
    values = CT_Colors
  ) +
  ggplot2::theme_classic(base_size = 15) +
  ggplot2::scale_y_log10(breaks = c(100, 1000, 10000, 100000), labels = c("100", "1k", "10k", "100k")) +
  ggplot2::geom_hline(yintercept = nCount_RNA_max, size = 1, linetype = 2) +
  ggplot2::geom_hline(yintercept = nCount_RNA_min, size = 1, linetype = 2) +
  ggplot2::theme(
    legend.position  = "",
    axis.title       = ggplot2::element_blank(),
    axis.text.x      = ggplot2::element_text(angle = 45, hjust = 1),
    strip.background = ggplot2::element_blank()
  )
ggplot2::ggsave(
  filename = paste0(plot_path, "DataQuality/Counts_CT.png"),
  width = 6,
  height = 5.3
)

# Genes per Cell
quality_metric <- "nFeature_RNA"
data <- tidyr::as_tibble(t(as.matrix(
  ds@meta.data[[quality_metric]]
)))
data <- data.frame("Value" = ds@meta.data[[quality_metric]])
data$col <- ds@meta.data$Celltype
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
  ggplot2::geom_violin(
    scale = "width", size = 0.5, alpha = 0.75, draw_quantiles = 0.5,
    trim = TRUE
  ) +
  ggplot2::scale_fill_manual(
    values = CT_Colors
  ) +
  ggplot2::theme_classic(base_size = 15) +
  ggplot2::scale_y_log10(breaks = c(500, 1000, 5000, 10000), labels = c("500", "1,000", "5,000", "10,000")) +
  ggplot2::geom_hline(yintercept = nFeature_RNA_max, linetype = 2, size = 1) +
  ggplot2::geom_hline(yintercept = nFeature_RNA_min, linetype = 2, size = 1) +
  ggplot2::theme(
    legend.position  = "",
    axis.title       = ggplot2::element_blank(),
    axis.text.x      = ggplot2::element_text(angle = 45, hjust = 1),
    strip.background = ggplot2::element_blank()
  )
ggplot2::ggsave(
  filename = paste0(plot_path, "DataQuality/Genes_CT.png"),
  width = 6,
  height = 5.3
)

# mitochondrial count fraction
quality_metric <- "percent.mt"
data <- tidyr::as_tibble(t(as.matrix(
  ds@meta.data[[quality_metric]]
)))
data <- data.frame("Value" = ds@meta.data[[quality_metric]])
data$col <- ds@meta.data$Celltype
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
  ggplot2::geom_violin(
    scale = "width", size = 0.5, alpha = 0.75, draw_quantiles = 0.5,
    trim = TRUE
  ) +
  ggplot2::scale_fill_manual(
    values = CT_Colors
  ) +
  ggplot2::theme_classic(base_size = 15) +
  ggplot2::geom_hline(yintercept = percent.mt_max, linetype = 2, size = 1) +
  ggplot2::geom_hline(yintercept = percent.mt_min, linetype = 2, size = 1) +
  ggplot2::theme(
    legend.position  = "",
    axis.title       = ggplot2::element_blank(),
    axis.text.x      = ggplot2::element_text(angle = 45, hjust = 1),
    strip.background = ggplot2::element_blank()
  )
ggplot2::ggsave(
  filename = paste0(plot_path, "DataQuality/Mito_CT.png"),
  width = 6,
  height = 5.3
)


# Make Plot for CD4/CD8 ratio
data <- data.frame(
  table(ds@meta.data$Celltype[ds@meta.data$Celltype %in% c("T CD4", "T CD8", "Treg")])
)
colnames(data) <- c("Annotation", "Abundance")
data <- data[c(8,9,11),]
ratio <- round(sum(data$Abundance[data$Annotation == "Treg"],
                   data$Abundance[data$Annotation == "T CD4"])/data$Abundance[data$Annotation == "T CD8"], 2)

# Compute the position of labels
data <- data %>% 
  arrange(desc(Annotation)) %>%
  mutate(prop = Abundance / sum(data$Abundance) *100) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )

# Basic piechart
plot <- ggplot2::ggplot(data, ggplot2::aes(x="", y=prop, fill=Annotation)) +
  ggplot2::geom_bar(stat="identity", width=1, color="white") +
  ggplot2::coord_polar("y", start=0) +
  ggplot2::theme_void() +
  ggplot2::theme(legend.position="none") +
  ggplot2::geom_text(ggplot2::aes(y = ypos, label = Abundance), color = "black", size=7) +
  ggplot2::scale_color_manual(values = CT_Colors) +
  ggplot2::scale_fill_manual(values = CT_Colors) 

plot <- plot + 
  ggplot2::geom_text(ggplot2::aes(x = -0.2, label = paste0("CD4/CD8-ratio\n = \n", ratio)), color = "black", size = 7)

ggplot2::ggsave(
  plot = plot,
  filename = paste0(plot_path, "DataQuality/CD4CD8ratio.png"),
  width = 7,
  height = 7
)

saveRDS(ds, file = "/home/alexander/Data/members/Alexander_Leipold/projects/Leo_Rasche/MM_Case_Report/Revision/13_05_22/CP_BAL_Annotated.Rds", compress = FALSE)

###############################################################################
# Extract T cells
ds <- subset(ds, subset = Celltype == "T CD4" | Celltype == "T CD8" |
                    Celltype == "Treg")

# Normalize and find highly variable genes
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
ds <- Seurat::RunPCA(ds, npcs = 55)
ds <- Seurat::RunUMAP(
  object    = ds,
  dims      = 1:15,
  reduction = "pca",
  seed.use  = seed
)

# Clustering with Louvain algorithm based on SNN-graph
ds <- Seurat::FindNeighbors(
  object    = ds,
  dims      = 1:15,
  reduction = "pca"
)
ds <- Seurat::FindClusters(
  object     = ds,
  resolution = 0.7
)

# Plot clusters
embedding <- "umap"
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
    size = 0.9
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
  filename = paste0(plot_path, "CP_BAL_T_cells/Cluster_CP_BAL_T.png"),
  width = 7,
  height = 5
)

# Plot Quality per clusters
# UMI Counts per cell
quality_metric <- "nCount_RNA"
data <- tidyr::as_tibble(t(as.matrix(
  ds@meta.data[[quality_metric]]
)))
data <- data.frame("Value" = ds@meta.data[[quality_metric]])
data$col <- ds@meta.data$seurat_clusters
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
  ggplot2::geom_violin(
    scale = "width", size = 0.5, alpha = 0.75, draw_quantiles = 0.5,
    trim = TRUE
  ) +
  ggplot2::theme_classic(base_size = 15) +
  ggplot2::scale_y_log10(breaks = c(2000, 5000, 10000, 20000), labels = c("2k", "5k", "10k", "20k")) +
  ggplot2::theme(
    legend.position  = "",
    axis.title       = ggplot2::element_blank(),
    axis.text.x      = ggplot2::element_text(angle = 0),
    strip.background = ggplot2::element_blank()
  ) 
ggplot2::ggsave(
  filename = paste0(plot_path, "CP_BAL_T_cells/QC_counts.png"),
  width = 6,
  height = 5.3
)

# Genes per Cell
quality_metric <- "nFeature_RNA"
data <- tidyr::as_tibble(t(as.matrix(
  ds@meta.data[[quality_metric]]
)))
data <- data.frame("Value" = ds@meta.data[[quality_metric]])
data$col <- ds@meta.data$seurat_clusters
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
  ggplot2::geom_violin(
    scale = "width", size = 0.5, alpha = 0.75, draw_quantiles = 0.5,
    trim = TRUE
  ) +
  ggplot2::theme_classic(base_size = 15) +
  ggplot2::scale_y_log10(breaks = c(500, 1000, 2000, 3000, 5000, 10000), labels = c("500", "1,000", "2,000", "3,000", "5,000", "10,000")) +
  ggplot2::theme(
    legend.position  = "",
    axis.title       = ggplot2::element_blank(),
    axis.text.x      = ggplot2::element_text(angle = 0),
    strip.background = ggplot2::element_blank()
  )
ggplot2::ggsave(
  filename = paste0(plot_path, "CP_BAL_T_cells/QC_genes.png"),
  width = 6,
  height = 5.3
)

# mitochondrial count fraction
quality_metric <- "percent.mt"
data <- tidyr::as_tibble(t(as.matrix(
  ds@meta.data[[quality_metric]]
)))
data <- data.frame("Value" = ds@meta.data[[quality_metric]])
data$col <- ds@meta.data$seurat_clusters
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
  ggplot2::geom_violin(
    scale = "width", size = 0.5, alpha = 0.75, draw_quantiles = 0.5,
    trim = TRUE
  ) +
  ggplot2::theme_classic(base_size = 15) +
  ggplot2::theme(
    legend.position  = "",
    axis.title       = ggplot2::element_blank(),
    axis.text.x      = ggplot2::element_text(angle = 0),
    strip.background = ggplot2::element_blank()
  )
ggplot2::ggsave(
  filename = paste0(plot_path, "CP_BAL_T_cells/QC_mito.png"),
  width = 6,
  height = 5.3
)

# DE analysis Clusters
ds <- Seurat::SetIdent(ds, value = ds@meta.data$seurat_clusters)
markers <- Seurat::FindAllMarkers(ds, only.pos = TRUE)
markers$gene <- to_genes[markers$gene]
markers <- markers[,c(6, 7, 2, 3, 4, 1, 5)]
markers$ENSEMBL_ID <- to_ids[markers$gene]
# Save as csv
write.csv(markers, file = "/home/alexander/Data/members/Alexander_Leipold/projects/Leo_Rasche/MM_Case_Report/Revision/13_05_22/tables/CP_BAL_T_Cluster_DE.csv", row.names = FALSE)


# Plot marker/important genes
pltlist <- list()
genes <- c(
  "IFNG", "TNF", "RORC", "TBX21", "CCR6", "CXCR3", "KLRB1",
  "IL23R", "IL4I1", "ABCB1", "CXCR6", "CSF2", "CCL20", "IL17A",
  
  "CD4", "CD40LG", "CD8A", "CD8B", "STAT3", "RUNX1", "RUNX3", "IL17B", "IL25", "IL17F",
  "IL26", "IL2", "CCR7", "TCF7", "S1PR1", "NKG7", "GZMK", "PRF1", "FOXP3", "IL2RA"
)

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
      size = 0.85
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

pltlistFig2 <- list(plot_IFNG, plot_TNF, plot_CSF2, plot_RORC, plot_TBX21, plot_CCR6, plot_CXCR3,
                    voidplt, voidplt, voidplt, voidplt, voidplt, voidplt, voidplt,
                    plot_KLRB1, plot_IL23R, plot_IL4I1, plot_ABCB1, plot_CXCR6, plot_CCL20, plot_IL17A)
plot <- ggpubr::ggarrange(plotlist = pltlistFig2,
                          nrow = 3, ncol = 7,
                          heights = c(1, -0.6, 1))
ggplot2::ggsave(
  plot = plot,
  filename = paste0(plot_path, "CP_BAL_T_cells/CP_BAL_T_GenesA.png"),
  width = 16,
  height = 8.5
)

# UMAPxgene plots Figure S3
pltlistFigS3 <- list(plot_CD4, plot_CD40LG, plot_CD8A, plot_CD8B, plot_STAT3, plot_RUNX1, plot_RUNX3, plot_IL17B, plot_IL25, plot_IL17F,
                     voidplt, voidplt, voidplt, voidplt, voidplt, voidplt, voidplt, voidplt, voidplt, voidplt,
                    plot_IL26, plot_IL2, plot_CCR7, plot_TCF7, plot_S1PR1, plot_NKG7, plot_GZMK, plot_PRF1, plot_FOXP3, plot_IL2RA)
plot <- ggpubr::ggarrange(plotlist = pltlistFigS3,
                          nrow = 3, ncol = 10,
                          heights = c(1, -0.6, 1))

ggplot2::ggsave(
  plot = plot,
  filename = paste0(plot_path, "CP_BAL_T_cells/CP_BAL_T_GenesB.png"),
  width = 23,
  height = 8.5
)


# Heatmap of genes for T subset annotation--------------------------
data <- data.frame("Cluster" = ds@meta.data$seurat_clusters,
                   "CD4" = ds@assays$RNA@data[to_ids["CD4"],],
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
                   "CD40LG" = ds@assays$RNA@data[to_ids["CD40LG"],],
                   "IL7R" = ds@assays$RNA@data[to_ids["IL7R"],],
                   "SELL" = ds@assays$RNA@data[to_ids["SELL"],],
                   "NKG7" = ds@assays$RNA@data[to_ids["NKG7"],],
                   "PRF1" = ds@assays$RNA@data[to_ids["PRF1"],],
                   "FOXP3" = ds@assays$RNA@data[to_ids["FOXP3"],],
                   "IL2RA" = ds@assays$RNA@data[to_ids["IL2RA"],],
                   "CCR7" = ds@assays$RNA@data[to_ids["CCR7"],],
                   "TCF7" = ds@assays$RNA@data[to_ids["TCF7"],],
                   "S1PR1" = ds@assays$RNA@data[to_ids["S1PR1"],]
                   )
library(dplyr)
data <- data %>% group_by(Cluster) %>% summarise(across(everything(), list(mean)))
colnames(data) <- c("Cluster", gsub("_1", "", colnames(data[, 2:length(colnames(data))])))
data <- as.matrix(data[, 2:length(colnames(data))])
data <- scale(data)
data[data > 1.5] <- 1.5
data[data < -1] <- -1

map <- pheatmap::pheatmap(mat = data,
                          color = colorRampPalette(rev(RColorBrewer::brewer.pal(n=7, name = "RdBu")))(500)[100:500],
                          scale = "none",
                          cluster_rows = F,
                          cluster_cols = F
                          )
ggplot2::ggsave(
  plot = map,
  filename = paste0(plot_path, "CP_BAL_T_cells/CP_BAL_T_Heatmap_v2.png"),
  width = 8.3,
  height = 3
)

# test genuinity per cell by Coexpression analysis----------
#data <- data.frame("Cell" = rownames(ds@meta.data),
#                   "Cluster" = ds@meta.data$seurat_clusters,
#                   "Subset" = ds@meta.data$T_subset,
#                   "TBX21" = ds@assays$RNA@data[to_ids["TBX21"],],
#                   "IFNG" = ds@assays$RNA@data[to_ids["IFNG"],],
#                   "CXCR3" = ds@assays$RNA@data[to_ids["CXCR3"],],
#                   "RORC" = ds@assays$RNA@data[to_ids["RORC"],],
#                   "CCL20" = ds@assays$RNA@data[to_ids["CCL20"],],
#                   "CCR6" = ds@assays$RNA@data[to_ids["CCR6"],])
# get cells that are coexpressing one of the gene-pairs
# RORC/TBX21 -> transcription factors - TF
# CCL20/IFNG -> cytokines - Cyt
# CCR6/CXCR3 -> Receptors - rec


data <- data.frame("Cell" = rownames(ds@meta.data),
                   "Cluster" = ds@meta.data$seurat_clusters,
                   "Subset" = ds@meta.data$T_subset,
                   "TBX21" = ds@assays$RNA@data[to_ids["TBX21"],],
                   "IFNG" = ds@assays$RNA@data[to_ids["IFNG"],],
                   "CXCR3" = ds@assays$RNA@data[to_ids["CXCR3"],],
                   "TNF" = ds@assays$RNA@data[to_ids["TNF"],],
                   "RORC" = ds@assays$RNA@data[to_ids["RORC"],],
                   "CCL20" = ds@assays$RNA@data[to_ids["CCL20"],],
                   "CCR6" = ds@assays$RNA@data[to_ids["CCR6"],],
                   "IL23R" = ds@assays$RNA@data[to_ids["IL23R"],],
                   "KLRB1" = ds@assays$RNA@data[to_ids["KLRB1"],])
data$Coexpr <- "No Co-expression"
data$Coexpr[data$TBX21 + data$IFNG + data$CXCR3 + data$TNF > 0 & data$RORC + data$CCL20 + data$CCR6 + data$IL23R > 0] <- "Co-expression"

ds@meta.data[data$Cell,]$Coexprez <- data$Coexpr

# save csv
write.csv(markers, file = "/home/alexander/Data/members/Alexander_Leipold/projects/Leo_Rasche/MM_Case_Report/Revision/13_05_22/tables/CP_BAL_T_Coexpr.csv", row.names = FALSE)


#double_pos_TF_cells <- as.character(data$Cell[data$RORC > 0 & data$TBX21 > 0])
#double_pos_Cytokine_cells <- as.character(data$Cell[data$CCL20 > 0 & data$IFNG > 0])
#double_pos_Receptor_cells <- as.character(data$Cell[data$CCR6 > 0 & data$CXCR3 > 0])

# genuine by coexpression; exact
# add meta.data colums for coexpression of gene pairs - boolean
#ds@meta.data$genuineRec <- "False"
#ds@meta.data$genuineRec[rownames(ds@meta.data) %in% double_pos_Receptor_cells] <- "True"
#ds@meta.data$genuineCyt <- "False"
#ds@meta.data$genuineCyt[rownames(ds@meta.data) %in% double_pos_Cytokine_cells] <- "True"
#ds@meta.data$genuineTF <- "False"
#ds@meta.data$genuineTF[rownames(ds@meta.data) %in% double_pos_TF_cells] <- "True"##

#ds@meta.data$Coexpr <- "No Co-expression"
#ds@meta.data$Coexpr[ds@meta.data$genuineTF == "True"] <- "Single TF"
#ds@meta.data$Coexpr[ds@meta.data$genuineRec == "True"] <- "Single Rec"
#ds@meta.data$Coexpr[ds@meta.data$genuineCyt == "True"] <- "Single Cyt"
#ds@meta.data$Coexpr[ds@meta.data$genuineCyt == "True" &
#                           ds@meta.data$genuineTF == "True"] <- "Double-TF+Cyt"
#ds@meta.data$Coexpr[ds@meta.data$genuineRec == "True" &
#                           ds@meta.data$genuineTF == "True"] <- "Double-TF+Rec"
#ds@meta.data$Coexpr[ds@meta.data$genuineRec == "True" &
#                           ds@meta.data$genuineCyt == "True"] <- "Double-Cyt+Rec"
#ds@meta.data$Coexpr[ds@meta.data$genuineRec == "True" &
#                           ds@meta.data$genuineCyt == "True" &
#                           ds@meta.data$genuineTF == "True"] <- "Triple"
#ds@meta.data$Coexpr <- as.factor(ds@meta.data$Coexpr)
#ds@meta.data$Coexpr <- factor(ds@meta.data$Coexpr, levels(ds@meta.data$Coexpr)[c(4, 6, 5, 7, 1, 3, 2, 8)])

#ds@meta.data$Coexprez <- "Co-expression"
#ds@meta.data$Coexprez[ds@meta.data$Coexpr == "No Co-expression"] <- "No Co-expression"

td <- table(ds@meta.data$seurat_clusters, ds@meta.data$Coexprez)
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
data$CoexprType <- factor(data$CoexprType, levels(data$CoexprType)[c(2, 1)])
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
  ggplot2::scale_color_manual(
    values = c("gray65", "red")
  ) +
  ggplot2::scale_fill_manual(
    values = c("gray65", "red")
  )
ggplot2::ggsave(
  plot = plot,
  filename = paste0(plot_path, "CP_BAL_T_cells/CP_BAL_T_Coexprez_Cluster.png"),
  width = 6.5,
  height = 5.5
)

# Differential expression between clusters
markers <- Seurat::FindAllMarkers(ds, only.pos = T )
markers$gene <- to_genes[markers$gene]
View(markers)

# T-subset annotation
cluster.annotation <- c(
  "0" = "Th17.1",
  "1" = "Th17.1",
  "2" = "CD4 low-count/gene",
  "3" = "CD8 Th17.1-like",
  "4" = "CD8 CTL",
  "5" = "Treg",
  "6" = "CD4 TCM-like"
)
ds@meta.data$T_subset <- factor(cluster.annotation[ds@meta.data$seurat_clusters], 
                                     unique(cluster.annotation)[c(1, 3, 2, 6, 4, 5)])


# DE analysis T_subset
ds <- Seurat::SetIdent(ds, value = ds@meta.data$T_subset)
markers <- Seurat::FindAllMarkers(ds, only.pos = TRUE)
markers$gene <- to_genes[markers$gene]
markers <- markers[,c(6, 7, 2, 3, 4, 1, 5)]
markers$ENSEMBL_ID <- to_ids[markers$gene]
# Save as csv
write.csv(markers, file = "/home/alexander/Data/members/Alexander_Leipold/projects/Leo_Rasche/MM_Case_Report/Revision/13_05_22/tables/CP_BAL_T_Tsubset_DE.csv", row.names = FALSE)
ds <- Seurat::SetIdent(ds, value = ds@meta.data$seurat_clusters)


# Subset colors
subset.colors <- c(
  "Th17.1" = "skyblue",
  "CD8 Th17.1-like" = "steelblue3",
  "CD4 low-count/gene" = "darkseagreen",
  "CD4 TCM-like" = "palegreen",
  "CD8 CTL" = "lightsalmon2",
  "Treg" = "skyblue4"
)

# Plot Subsets
embedding <- "umap"
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
    size = 0.9
  ) +
  ggplot2::scale_color_manual(values = subset.colors) +
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
  filename = paste0(plot_path, "CP_BAL_T_cells/Subset_CP_BAL_T.png"),
  width = 6.5,
  height = 5
)

# Plot UMAP with genuinity
embedding <- "umap"
data <- tidyr::as_tibble(
  x = ds@reductions[[embedding]]@cell.embeddings
)
colnames(data) <- c("x", "y")
category <- "Coexprez"
data$col <- as.factor(ds@meta.data[[category]])
data$col <- factor(data$col, levels(data$col)[c(2,1)])
plot <- ggplot2::ggplot(
  data    = data[data$col == "Co-expression", ],
  mapping = ggplot2::aes(
    x   = x,
    y   = y,
    col = col
  )
) +
  ggplot2::geom_point(
    data = data[data$col == "No Co-expression", ],
    ggplot2::aes(x, y),
    color = "gray65", size = 0.85) +
  ggplot2::geom_point(
    size = 1,
    color = "red"
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
  filename = paste0(plot_path, "CP_BAL_T_cells/genuinity_UMAP.png"),
  width = 5.5,
  height = 5
)

# Plot genuinity per T-subset
td <- table(ds@meta.data$T_subset, ds@meta.data$Coexprez)
tdpct <- td_to_pct(td)

data <- data.frame(tdpct)
colnames(data) <- c("Subset", "CoexprType", "proportion")
data$CoexprType <- factor(data$CoexprType, levels(data$CoexprType)[c(2, 1)])
plot <- plot <- ggplot2::ggplot(
  data    = data,
  mapping = ggplot2::aes(
    x     = Subset,
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
    axis.text.x      = ggplot2::element_text(angle = 45, hjust = 1), 
    axis.title       = ggplot2::element_blank(),
    legend.position  = "right",
    legend.text      = ggplot2::element_text(size = 20),
    legend.title     = ggplot2::element_blank()
  ) +
  ggplot2::scale_color_manual(values = c("gray65", "red")) +
  ggplot2::scale_fill_manual(values = c("gray65", "red"))
ggplot2::ggsave(
  plot = plot,
  filename = paste0(plot_path, "CP_BAL_T_cells/CP_BAL_T_Coexprez_Subset.png"),
  width = 6,
  height = 6
)


# Violin plots of IFNG/TNF/RORC/CD40LG
# Genes per Cell

assay <- "RNA"
genes <- c("IFNG", "TNF", "KLRB1", "CD40LG")
seed <- 1999
data <- tidyr::as_tibble(t(as.matrix(
  ds@assays[[assay]]@data[to_ids[genes],]
)))
colnames(data) <- genes
data$col <- as.factor(ds@meta.data$T_subset)
data1 <- data
for (i in genes) {
  noise <- rnorm(n = nrow(data[,i])) / 100000
  data[,i] <- data[,i] + noise 
}
data <- tidyr::gather(data, "Gene", "Value", -col)
data$Gene <- factor(data$Gene, unique(data$Gene))
plot <- ggplot2::ggplot(
  data    = data,
  mapping = ggplot2::aes(
    x   = col,
    y   = Value,
    fill = col
  )
) +
  ggplot2::geom_point(
    position = "jitter", size = 0.3, color = "black", 
  ) +
  ggplot2::geom_violin(
    scale = "width", size = 1, alpha = 0.75, draw_quantiles = 0.5,
    trim = TRUE
  ) +
  ggplot2::scale_fill_manual(values = subset.colors) +
  ggplot2::scale_color_distiller(palette = "Reds", direction = 1) +
  ggplot2::theme_classic(base_size = 15) +
  ggplot2::facet_wrap(~Gene, scales = "free_y", nrow = 1) +
  ggplot2::theme(
    legend.position  = "",
    axis.title       = ggplot2::element_blank(),
    axis.text.x      = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
    strip.text.x     = ggplot2::element_text(face = "italic", size = 25),
    strip.text.y     = ggplot2::element_text(angle = 0, face = "plain"),
    strip.background = ggplot2::element_blank(),
    panel.border     = ggplot2::element_rect(
      fill = NA, colour = "black", size = 0.3
    )
  )
data1 <- tidyr::gather(data1, "Gene", "Value", -col,)
data1$Gene <- factor(data1$Gene, unique(data1$Gene))
my_comparisons <- list(c("Th17.1", "CD8 Th17.1-like"))
plot <- plot + ggpubr::stat_compare_means(
  data = data1,
  comparisons = my_comparisons,
  method = "wilcox", label = "p.format"
) + ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, .1)))
ggplot2::ggsave(
  plot = plot,
  filename = paste0(plot_path, "CP_BAL_T_cells/CP_BAL_T_GeneViolins.png"),
  width = 15,
  height = 7.5
)

saveRDS(ds, file = "/home/alexander/Data/members/Alexander_Leipold/projects/Leo_Rasche/MM_Case_Report/Revision/13_05_22/CP_BAL_T_Annotated.Rds", compress = FALSE)






