# Script for analysis of EMD sample
# Alexander Leipold


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
  filename = paste0(plot_path, "EMD/CellType.png"),
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
  filename = paste0(plot_path, "EMD/GenesA.png"),
  width = 10,
  height = 10
)












#saveRDS(ds, file = "/home/alexander/Data/members/Alexander_Leipold/projects/Leo_Rasche/MM_Case_Report/Revision/13_05_22/EMD.Rds")




plot <- Seurat::DimPlot(ds, group.by = "Celltype")
plot

Seurat::FeaturePlot(ds, features = to_ids["HIF1A"], order = T, pt.size = 1)
Seurat::FeaturePlot(ds, features = to_ids["CD3D"], order = T, pt.size = 1)



ds@meta.data[selecia,]$Celltype <- "Fibroblast"


ds_MMs <- subset(ds_MMs, subset = Type != "EMD")

ds <- merge(ds, y = ds_MMs)

########################################################################

# Script for comparison of EMD and IMD MM transcriptomes
# Alexander Leipold

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




colorsPH <- RColorBrewer::brewer.pal(6, "Greens")
colorsPT <- RColorBrewer::brewer.pal(7, "Purples")
Origin_Colors <- c("EMD" = "darkorange2",
                   "DaVia_Relapse"  = "orchid",
                   "DeJong_PH1" = colorsPH[1],
                   "DeJong_PH2" = colorsPH[2],
                   "DeJong_PH3" = colorsPH[3],
                   "DeJong_PH4" = colorsPH[4],
                   "DeJong_PH5" = colorsPH[5],
                   "DeJong_PH6" = colorsPH[6],
                   "DeJong_PT1" = colorsPT[1],
                   "DeJong_PT2" = colorsPT[2],
                   "DeJong_PT3" = colorsPT[3],
                   "DeJong_PT4" = colorsPT[4],
                   "DeJong_PT5" = colorsPT[5],
                   "DeJong_PT6" = colorsPT[6],
                   "DeJong_PT7" = colorsPT[7])

Seurat::DimPlot(ds, group.by = "Sample", cols = Origin_Colors)


Type_Colors <- c("EMD" = "darkorange2",
                 "IMD \npost CAR T" = "orchid",
                 "IMD \nnewly diagnosed" = "aquamarine2")

Seurat::DimPlot(ds, group.by = "Type", cols = Type_Colors)

############################
embedding <- "umap"
data <- tidyr::as_tibble(
  x = ds@reductions[[embedding]]@cell.embeddings
)
colnames(data) <- c("x", "y")
category <- "Sample"
data$col <- as.factor(ds@meta.data[[category]])
data$col <- factor(data$col, levels(data$col)[c(15,seq(1,14))])
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
      override.aes = list(size = 5),
      reverse = TRUE
    )) +
  ggplot2::scale_color_manual(
    values = Origin_Colors[rev(levels(data$col))]
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

plot

filename <- paste0(
  plot_path,
  "EMD_IMD/",
  embedding, "_", category,
  ".png"
)

ggplot2::ggsave(
  filename = filename,
  plot     = plot,
  width    = 9.6,
  height   = 6
)

# 6.2 UMAP Type----------------------------------------------------------------
category <- "Type"
data$col <- as.factor(ds@meta.data[[category]])
data$col <- factor(data$col, levels(data$col)[c(1, 3, 2)])
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
      override.aes = list(size = 5),
      reverse = TRUE
    )) +
  ggplot2::scale_color_manual(
    values = Type_Colors[rev(levels(data$col))]
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

plot

filename <- paste0(
  plot_path,
  "EMD_IMD/",
  embedding, "_", category,
  ".png"
)

ggplot2::ggsave(
  filename = filename,
  plot     = plot,
  width    = 9.65,
  height   = 6
)

##############
assay <- "RNA"
genes <- c("EPCAM", "SFN", "KRT8", "KRT18", 
           "ADGRE2", "ITGB7", "VEGFA", "PTN",
           "ANG", "TGFA")
seed <- 1999
data <- tidyr::as_tibble(t(as.matrix(
  ds@assays[[assay]]@data[to_ids[genes],]
)))
colnames(data) <- genes
data$Condition <- as.factor(ds@meta.data$Type)
#data$Origin <- as.factor(ds@meta.data$Origin)
data$Condition <- factor(data$Condition, levels(data$Condition)[c(1, 3, 2)])
#rndsample1 <- dplyr::sample_n(data[data$Condition == "EMD", ], size = 1300)
#rndsample2 <- dplyr::sample_n(data[data$Condition == "IMD \npost CAR T", ], size = 1300)
#rndsample3 <- dplyr::sample_n(data[data$Condition == "IMD \nnewly diagnosed", ],size = 1000)
#rndsample3 <- dplyr::sample_n(data[data$Origin == "PH1", ],size = 100)
#rndsample4 <- dplyr::sample_n(data[data$Origin == "PH2", ],size = 100)
#rndsample5 <- dplyr::sample_n(data[data$Origin == "PH3", ],size = 100)
#rndsample6 <- dplyr::sample_n(data[data$Origin == "PH4", ],size = 100)
#rndsample7 <- dplyr::sample_n(data[data$Origin == "PH5", ],size = 100)
#rndsample8 <- dplyr::sample_n(data[data$Origin == "PH6", ],size = 100)
#rndsample9 <- dplyr::sample_n(data[data$Origin == "PT1", ],size = 100)
#rndsample10 <- dplyr::sample_n(data[data$Origin == "PT2_own", ],size = 100)
#rndsample11 <- dplyr::sample_n(data[data$Origin == "PT3", ],size = 100)
#rndsample12 <- dplyr::sample_n(data[data$Origin == "PT4", ],size = 100)
#rndsample13 <- dplyr::sample_n(data[data$Origin == "PT5", ],size = 100)
#rndsample14 <- dplyr::sample_n(data[data$Origin == "PT6", ],size = 100)
#rndsample15 <- dplyr::sample_n(data[data$Origin == "PT7", ],size = 100)

#data <- rbind(rndsample1, rndsample2, rndsample3)
#data <- rbind(rndsample1, rndsample2, rndsample3, rndsample4, rndsample5, rndsample6,
#              rndsample7, rndsample8, rndsample9, rndsample10, rndsample11,
#              rndsample12, rndsample13, rndsample14, rndsample15)
#data$Origin <- NULL


data1 <- data
for (i in genes) {
  noise <- rnorm(n = nrow(data[,i])) / 100000
  data[,i] <- data[,i] + noise 
}

#data$Condition <- factor(data$Condition, levels(data$Condition)[c(1,3,2)])
data <- tidyr::gather(data, "Gene", "Value", -Condition)
data$Gene <- factor(data$Gene, unique(data$Gene))



plot <- ggplot2::ggplot(
  data    = data,
  mapping = ggplot2::aes(
    x   = Condition,
    y   = Value,
    fill = Condition
  )
) +
  ggplot2::geom_point(
    position = "jitter", size = 0.3, color = "black", 
  ) +
  ggplot2::geom_violin(
    scale = "width", size = 1, alpha = 0.75, draw_quantiles = 0.5,
    trim = TRUE
  ) +
  ggplot2::scale_fill_manual(values = Type_Colors[levels(data$Condition)]) +
  ggplot2::scale_color_distiller(palette = "Reds", direction = 1) +
  ggplot2::theme_classic(base_size = 20) +
  ggplot2::facet_wrap(~Gene, scales = "free_y", nrow = 2) +
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
  )#  +
#ggplot2::scale_y_continuous(limits = c(0, 6))

data1 <- tidyr::gather(data1, "Gene", "Value", -Condition,)
data1$Gene <- factor(data1$Gene, unique(data1$Gene))
my_comparisons <- list(c("EMD", "IMD \npost CAR T"), c("EMD", "IMD \nnewly diagnosed"))
plot <- plot + ggpubr::stat_compare_means(
  data = data1,
  comparisons = my_comparisons,
  method = "wilcox", label = "p.format"
) + ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, .1)))

plot

filename <- paste0(
  plot_path,
  "EMD_IMD/Violins",
  ".png"
)
ggplot2::ggsave(
  filename = filename,
  plot     = plot,
  width    = 12,
  height   = 8
)


#############################
pltlist <- list()
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
      size = 0.5
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
  
  #if (i == tail(genes, n = 1)) {
  #  plot <- plot +
  #    ggplot2::theme(
  #      legend.position = "right"
  #    )
  #}
  
  assign(paste0("plot_", i), plot)
  pltlist[[paste0("plot_", i)]] <- plot
}


# put them together
voidplt <- ggplot2::ggplot() + ggplot2::theme_void()
pltlist2 <- list(plot_EPCAM, plot_SFN, plot_KRT8, plot_KRT18, plot_ADGRE2,
                 voidplt, voidplt, voidplt,
                 voidplt, voidplt, plot_ITGB7, plot_VEGFA, plot_PTN, plot_ANG, plot_TGFA)
plot <- ggpubr::ggarrange(plotlist = pltlist2,
                          nrow = 3, ncol = 5,
                          heights = c(1, -0.335, 1))

filename <- paste0(
  plot_path,
  "EMD_IMD/UMAPs",
  ".png"
)

ggplot2::ggsave(
  filename = filename,
  plot     = plot,
  width    = 14,
  height   = 10
)

#####################################


genesrev <- rev(genes)
#genes <- c("")
data <- t(as.matrix(ds@assays$RNA@data[to_ids[genesrev], ]))
#data <- scale(data)
data <- tidyr::as_tibble(data)
colnames(data) <- genesrev
data$Celltype <- factor(ds@meta.data$Sample)
data$Celltype <- factor(data$Celltype, levels(data$Celltype)[c(15, 1, seq(2,14))])
data$Celltype <- factor(data$Celltype, levels(data$Celltype)[seq(15,1)])
data <- tidyr::gather(data, "Gene", "Expression", -Celltype)
data$Gene <- factor(
  x      = data$Gene,
  levels = rev(unique(data$Gene))
)

data$Pct <- data$Expression > 0
data$Freq <- rep(1, length(data$Pct))

data <- dplyr::group_by(data, Gene)
data <- dplyr::mutate(data, Scaled = scale(Expression)[, 1])

data <- dplyr::group_by(data, Celltype, Gene)
data <- dplyr::summarise(
  data,
  Mean   = mean(Expression),
  Scaled = mean(Scaled),
  Pct    = sum(Pct)/sum(Freq)*100
)
colorscale <- RColorBrewer::brewer.pal(n = 9, name = "RdBu")
#data$Scaled[data$Scaled > 0.81] <- 0.81
#data$Scaled[data$Scaled < -0.81] <- 0.81


plot <- ggplot2::ggplot(
  data = data,
  mapping = ggplot2::aes(
    y = Celltype,
    x = Gene,
    size = Pct,
    col = Scaled
  )
) +
  ggplot2::geom_point() +
  ggplot2::scale_color_gradient2(midpoint=0, low=colorscale[length(colorscale)], mid="white",
                                 high=colorscale[1],
  ) +#labels = c("min", "max"), breaks = c(-2, 2)) +
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
    legend.title     = ggplot2::element_blank(),
    legend.text      = ggplot2::element_text(size = 15)
  ) +
  ggplot2::guides(
    color = ggplot2::guide_colorbar(
      barheight = 10, barwidth = 0.6, frame.colour = "black", ticks = FALSE, order = 1,
      
    ),
    size  = ggplot2::guide_legend(
      label.position = "right"
    )
  ) +
  ggplot2::scale_size_area(max_size = 8)# +

plot

filename <- paste0(
  plot_path,
  "EMD_IMD/DP",
  ".png"
)

ggplot2::ggsave(
  filename = filename,
  plot     = plot,
  width    = 7,
  height   = 5
)

