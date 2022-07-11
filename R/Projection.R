# Reference mapping of cutaneous Sarc T-cells onto integrated BAL T-cell embedding

# identify mapping anchors using reference.reduction = "pca", as recommended
mapping.anchors <- Seurat::FindTransferAnchors(reference = ds_integrated_BAL_T,
                                               query = ds,
                                               dims = 1:10,
                                               reference.reduction = "pca",
                                               k.filter = 100)
# Map teh query cells onto the reference umap and perform label transfer of T-cell subsets
ds <- Seurat::MapQuery(anchorset = mapping.anchors,
                       reference = ds_integrated_BAL_T,
                       query = ds,
                       refdata = list(T_subset = "T_subset"),
                       reference.reduction = "pca",
                       reduction.model = "mnn_umap_model")

# label cells with prediction score below 0.4 as unassigned
ds@meta.data$filtered.prediction <- ds@meta.data$predicted.T_subset
ds@meta.data$filtered.prediction[ds@meta.data$predicted.T_subset.score < 0.4] <- "unassigned" 

# Calculate Signature Module Scores
gene_ids_Nathan <- to_ids[c("CCR6", "DPP4", "CTSH", "LGALS3", "BHLHE40", "PDE4D")]
gene_ids_Gideon <- to_ids[c("RORA", "RORC", "RBPJ", "BHLHE40", "FURIN", "COTL1", "CCL20", "CCR6",
                            "IL23R", "CXCR3", "IFNG", "TNF")]
CT <- list()
CT$Score_Nathan <- gene_ids_Nathan
CT$Score_Gideon <- gene_ids_Gideon
ds <- Seurat::AddModuleScore(
  object = ds, features = CT, name = names(CT), nbin = 24
)

# Density plot
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

data <- data.frame("cell" = c(rownames(ds_integrated_BAL_T@meta.data), rownames(ds@meta.data)),
                   "x" = c(ds_integrated_BAL_T@reductions$mnn_umap@cell.embeddings[,1], ds@reductions$ref.umap@cell.embeddings[,1]),
                   "y" = c(ds_integrated_BAL_T@reductions$mnn_umap@cell.embeddings[,2], ds@reductions$ref.umap@cell.embeddings[,2]),
                   "Condition" = c(ds_integrated_BAL_T@meta.data$Condition, rep("SkinSarc", 4733)),
                   "Prediction" = c(ds_integrated_BAL_T@meta.data$T_subset, ds@meta.data$filtered.prediction),
                   "Confidence" = c(rep(NA, 13357), ds@meta.data$predicted.T_subset.score)
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
  filename = paste0(plot_path, "density.png"),
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
  ggplot2::scale_color_manual(values = labtransfer.colors[levels(factor(ds@meta.data$filtered.prediction))]) +
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
  filename = paste0(plot_path, "labtransfer.png"),
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
  filename = paste0(plot_path, "labtransferConf.png"),
  width = 8.3,
  height = 5
)


# ViolinPlot Gideon Th17/1 Score ~Subsets with more than 30 cells assigned to it
data <- data.frame("Value" = ds@meta.data$exTh17Score_Nathan1)
data$col <- data$Value
limits <- c(-0.5, 1)
data$col[data$col > max(limits)] <- max(limits)
data$col[data$col < min(limits)] <- min(limits)
data$cluster <- ds@meta.data$filtered.prediction
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
  filename = paste0(plot_path, "Score_Nathan_Violin_Subsets.png"),
  width = 3,
  height = 3
)


# ViolinPlot Gideon Th17/1 Score ~Subsets with more than 30 cells assigned to it
data <- data.frame("Value" = ds@meta.data$exTh17Score_Gideon2)
data$col <- data$Value
limits <- c(-0.4, 0.6)
data$col[data$col > max(limits)] <- max(limits)
data$col[data$col < min(limits)] <- min(limits)
data$cluster <- ds@meta.data$filtered.prediction
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
  filename = paste0(plot_path, "Score_Gideon_Violin_Subsets.png"),
  width = 3,
  height = 3
)


#saveRDS(ds, "savepath/Damsky_T_cells_projected.Rds", compress = FALSE)
