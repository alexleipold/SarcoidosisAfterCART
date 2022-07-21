# Script to precompute a tree and perform multi UniFrac distance analysis using this precomputed tree
# based on the scUniFrac package and code from LeBlanc et al. https://doi.org/10.1016/j.ccell.2022.02.016

plot_path <- "/path/for/plots/"

clust.tree <- function(object, method = "pca", group = "seurat_clusters",
                       dim_red_use = "mnn", pcs_use = NULL, genes_use = NULL) {
  require(Seurat)
  require(ape)
  
  if(method == "pca") {
    data_use <- data.frame(row.names = c(1:pcs_use))
    
    for (i in levels(object@meta.data[[group]])) {
      pca_tmp <- object@reductions[[dim_red_use]]@cell.embeddings[ds@meta.data$seurat_clusters == i, 1:pcs_use]
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

tree <- clust.tree(ds, method = "pca", group = "seurat_clusters", pcs_use = 55)



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

reslts <- sc.unifrac.multi(object = ds, group_by = "Condition", clusts_name = "seurat_clusters", tree = tree, perm_iters = 2, n_cores = 1)
reslts$distance

# order reslts properly
reslts$distance <- reslts$distance[c(2, 4, 3, 1), c(2, 4, 3, 1)]
mat <- as.matrix(reslts$distance)

# plot Heatmap
map <- pheatmap::pheatmap(mat = mat,
                          cluster_rows = F,
                          cluster_cols = F,
                          color = colorRampPalette(c("tomato", "lightyellow", "steelblue4"))(94))
)
ggplot2::ggsave(
  plot = map,
  filename = paste0(plot_path, "UniFrac.png"),
  width = 4.7,
  height = 4
)
