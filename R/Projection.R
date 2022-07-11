# Reference mapping of Sarc T-cells onto integrated BAL T-cell embedding-------

# get UMAP woth model
# UMAP with MNN-corrected PCA
ds_T <- Seurat::RunUMAP(
  object         = ds_T,
  dims           = 1:55,
  seed.use       = seed,
  reduction      = "mnn",
  return.model = TRUE,
  reduction.name = "mnn_umap_model")

mapping.anchors <- Seurat::FindTransferAnchors(reference = ds_T, query = ds_D_T, dims = 1:10, reference.reduction = "pca", k.filter = 100)
predictions <- Seurat::TransferData(anchorset = mapping.anchors, refdata = ds_T@meta.data$T_subset)
ds_D_T <- Seurat::MapQuery(anchorset = mapping.anchors, reference = ds_T, query = ds_D_T,
                           refdata = list(T_subset = "T_subset", clusters = "seurat_clusters", old_clust = "CP_Ann"),
                           reference.reduction = "pca", reduction.model = "mnn_umap_model")

p1 <- Seurat::DimPlot(ds_T, reduction = "mnn_umap", group.by = "T_subset", pt.size = 0.25)
p2 <- Seurat::DimPlot(ds_D_T, reduction = "ref.umap", group.by = "filtered.prediction", pt.size = 0.25)
p1+p2

p1 <- Seurat::DimPlot(ds_T, reduction = "mnn_umap", group.by = "seurat_clusters")
p2 <- Seurat::DimPlot(ds_D_T, reduction = "ref.umap", group.by = "predicted.clusters")
p1+p2

p1 <- Seurat::DimPlot(ds_T, reduction = "mnn_umap", group.by = "CP_Ann")
p2 <- Seurat::DimPlot(ds_D_T, reduction = "ref.umap", group.by = "predicted.old_clust")
p1+p2

p1 <- Seurat::DimPlot(ds_T, reduction = "mnn_umap", group.by = "T_subset")
p2 <- Seurat::DimPlot(ds_D_T, reduction = "ref.umap", group.by = "RNA_snn_res.0.5")
p1+p2

Seurat::FeaturePlot(ds_D_T, reduction = "ref.umap", features = to_ids["BHLHE40"], pt.size = 1, order = T)
Seurat::FeaturePlot(ds_D_T, reduction = "ref.umap", features = "predicted.T_subset.score", pt.size = 1, order = T)

gene_ids_Nathan <- to_ids[c("CCR6", "DPP4", "CTSH", "LGALS3", "BHLHE40", "PDE4D")]
gene_ids_Gideon <- to_ids[c("RORA", "RORC", "RBPJ", "BHLHE40", "FURIN", "COTL1", "CCL20", "CCR6",
                            "IL23R", "CXCR3", "IFNG", "TNF")]

CT <- list()
CT$exTh17Score_Nathan <- gene_ids_Nathan
CT$exTh17Score_Gideon <- gene_ids_Gideon

ds_D_T <- Seurat::AddModuleScore(
  object = ds_D_T, features = CT, name = names(CT), nbin = 24
)

Seurat::VlnPlot(subset(ds_D_T, subset = filtered.prediction %in% c("CD4 TCM ribosomal", "CD8 CTL", "Naive-like", "Th17.1", "Treg", "unassigned")),
                features = "exTh17Score_Nathan1", group.by = "filtered.prediction") +
  ggplot2::geom_violin(draw_quantiles = 0.5, alpha = 0.8, scale = "width")

Seurat::VlnPlot(ds_D_T,
                features = "exTh17Score_Gideon2", group.by = "filtered.prediction") +
  ggplot2::geom_violin(draw_quantiles = 0.5, alpha = 0.8, scale = "width")

Seurat::DimPlot(ds_D_T,
                reduction = "ref.umap", group.by = "Origin")

Seurat::FeaturePlot(ds_D_T, reduction = "ref.umap", features = "exTh17Score_Gideon2", pt.size = 1, order = T, cols = c("blue", "lightyellow", "red"))

Seurat::VlnPlot(ds_D_T, features = "predicted.T_subset.score", group.by = "predicted.T_subset")

ds_D_T@meta.data$filtered.prediction <- ds_D_T@meta.data$predicted.T_subset
ds_D_T@meta.data$filtered.prediction[ds_D_T@meta.data$predicted.T_subset.score < 0.4] <- "unassigned" 
table(ds_D_T@meta.data$filtered.prediction)



# Density plot for every Condition-----------------------
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

data <- data.frame("cell" = c(rownames(ds_T@meta.data), rownames(ds_D_T@meta.data)),
                   "x" = c(ds_T@reductions$mnn_umap@cell.embeddings[,1], ds_D_T@reductions$ref.umap@cell.embeddings[,1]),
                   "y" = c(ds_T@reductions$mnn_umap@cell.embeddings[,2], ds_D_T@reductions$ref.umap@cell.embeddings[,2]),
                   "Condition" = c(ds_T@meta.data$Condition, rep("SkinSarc", 4733)),
                   "Prediction" = c(ds_T@meta.data$T_subset, ds_D_T@meta.data$filtered.prediction),
                   "Confidence" = c(rep(NA, 13357), ds_D_T@meta.data$predicted.T_subset.score)
)

data$density <- NA
data[data$Condition == "SkinSarc",]$density <- scales::rescale(get_density(data[data$Condition == "SkinSarc",]$x, data[data$Condition == "SkinSarc",]$y, n = 100))

plot <- ggplot2::ggplot(
  data    = data[data$Condition == "SkinSarc", ],
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
  filename = paste0(plot_path, "SarcProjection/density.png"),
  width = 8.3,
  height = 5
)

labtransfer.colors <- c(
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
  "CD8 IFN-responsive" = "palevioletred2",
  "unassigned" = "gray83"
)

plot <- ggplot2::ggplot(
  data    = data[data$Condition == "SkinSarc", ],
  mapping = ggplot2::aes(
    x   = x,
    y   = y,
    col = Prediction
  )
) +
  ggplot2::geom_point(
    data = data,
    ggplot2::aes(x, y),
    color = "gray60", size = 3) +
  ggplot2::geom_point(
    size = 2
  ) +
  ggplot2::scale_color_manual(values = labtransfer.colors[levels(factor(ds_D_T@meta.data$filtered.prediction))]) +
  ggplot2::coord_fixed() +
  ggplot2::theme(
    panel.background = ggplot2::element_blank(),
    axis.ticks       = ggplot2::element_blank(),
    axis.text        = ggplot2::element_blank(),
    axis.title       = ggplot2::element_text(size = 15, hjust = 0.1, vjust = 0.1),
    legend.text      = ggplot2::element_text(size = 15),
    legend.title     = ggplot2::element_text(face = "plain", size = 18),
    legend.position  = "right",
    legend.key       = ggplot2::element_rect(fill = "white", color = "white"),
    plot.title            = ggplot2::element_text(face = "italic", size = 17, hjust = 0.5)
  ) +
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
ggplot2::ggsave(
  plot = plot,
  filename = paste0(plot_path, "SarcProjection/labtransfer.png"),
  width = 12,
  height = 5
)


data <- data[order(data$Confidence),]
data$Confidence[data$Confidence < 0.4] <- 0.4

plot <- ggplot2::ggplot(
  data    = data[data$Condition == "SkinSarc", ],
  mapping = ggplot2::aes(
    x   = x,
    y   = y,
    col = Confidence
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
    legend.position  = "right",
    legend.key       = ggplot2::element_rect(fill = "white", color = "white"),
    plot.title            = ggplot2::element_text(face = "italic", size = 17, hjust = 0.5)
  ) +
  ggplot2::guides(
    color = ggplot2::guide_legend(
      override.aes = list(size = 5),
      reverse = TRUE
    )) +
  ggplot2::scale_color_gradient2(
    low = "steelblue4",
    mid = "lightyellow",
    high = "red",
    midpoint = 0.7) +
  ggplot2::guides(
    color = ggplot2::guide_colorbar(
      barwidth  = 0.6,
      barheight = 10,
      ticks     = TRUE, 
      frame.colour = "black",
      title = "",
      label = TRUE
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
plot
ggplot2::ggsave(
  plot = plot,
  filename = paste0(plot_path, "SarcProjection/labtransferConf.png"),
  width = 8.3,
  height = 5
)


# ViolinPlot Gideon Th17/1 Score ~Subsets
data <- data.frame("Value" = ds_D_T@meta.data$exTh17Score_Nathan1)
data$col <- data$Value
limits <- c(-0.5, 1)
data$col[data$col > max(limits)] <- max(limits)
data$col[data$col < min(limits)] <- min(limits)
data$cluster <- ds_D_T@meta.data$filtered.prediction
data <- data[data$cluster %in% c("CD4 TCM ribosomal", "CD8 CTL", "Naive-like", "Th17.1", "Treg", "unassigned"), ]
data$cluster <- as.factor(data$cluster)
data$cluster <- factor(data$cluster, levels(data$cluster)[c(4, 1, 3, 5, 2, 6)])
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
  ggplot2::scale_fill_manual(values = labtransfer.colors) +
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
plot
ggplot2::ggsave(
  plot = plot,
  filename = paste0(plot_path, "SarcProjection/Score_Nathan_Violin_Subsets.png"),
  width = 3,
  height = 3
)


# ViolinPlot Gideon Th17/1 Score ~Subsets
data <- data.frame("Value" = ds_D_T@meta.data$exTh17Score_Gideon2)
data$col <- data$Value
limits <- c(-0.4, 0.6)
data$col[data$col > max(limits)] <- max(limits)
data$col[data$col < min(limits)] <- min(limits)
data$cluster <- ds_D_T@meta.data$filtered.prediction
data <- data[data$cluster %in% c("CD4 TCM ribosomal", "CD8 CTL", "Naive-like", "Th17.1", "Treg", "unassigned"), ]
data$cluster <- as.factor(data$cluster)
data$cluster <- factor(data$cluster, levels(data$cluster)[c(4, 1, 3, 5, 2, 6)])
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
  ggplot2::scale_fill_manual(values = labtransfer.colors) +
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
plot
ggplot2::ggsave(
  plot = plot,
  filename = paste0(plot_path, "SarcProjection/Score_Gideon_Violin_Subsets.png"),
  width = 3,
  height = 3
)






#saveRDS(ds_D_T, "/home/alexander/Data/members/Alexander_Leipold/projects/Leo_Rasche/MM_Case_Report/Revision/13_05_22/Damsky_T_cells_projected.Rds", compress = FALSE)
