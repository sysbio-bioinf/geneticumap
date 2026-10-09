options(Seurat.warn.umap.uwot = FALSE)

sc_clustering_pipeline <- function(count_data, config, params, return_markers = FALSE) {
  #check data is named list or matrix
  checkmate::assert(checkmate::check_list(count_data, names = "unique"), checkmate::check_class(count_data, "dgCMatrix"), combine = "or", .var.name = "count_data")

  seurat_obj <- Seurat::CreateSeuratObject(counts = count_data, project = "geneticumap", min.cells = 3, min.features = 200)

  # quality control and cell selection
  seurat_obj[["percent_mt"]] <- Seurat::PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > config$sc_clustering$quality_thresholds$nFeature_RNA_min &
    nFeature_RNA < config$sc_clustering$quality_thresholds$nFeature_RNA_max & percent_mt < config$sc_clustering$quality_thresholds$percent_mt)

  # round integer-valued parameters (optimizer works with continuous floats)
  params$nfeatures <- round(params$nfeatures)
  params$dimensionality <- round(params$dimensionality)
  params$umap_neighbors <- round(params$umap_neighbors)

  # normalization
  seurat_obj <- Seurat::NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = params$scale_factor, verbose = FALSE)

  # identification of highly variable features
  seurat_obj <- Seurat::FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = params$nfeatures, verbose = FALSE)

  # linear transformation (scaling)
  all.genes <- rownames(seurat_obj)
  seurat_obj <- Seurat::ScaleData(seurat_obj, features = all.genes, verbose = FALSE)

  # PCA
  seurat_obj <- Seurat::RunPCA(seurat_obj, features = Seurat::VariableFeatures(object = seurat_obj), verbose = FALSE)

  # cluster the cells
  seurat_obj <- Seurat::FindNeighbors(seurat_obj, dims = seq_len(params$dimensionality), verbose = FALSE)
  seurat_obj <- Seurat::FindClusters(seurat_obj, resolution = params$resolution, verbose = FALSE)

  # run UMAP
  seurat_obj <- Seurat::RunUMAP(seurat_obj, dims = seq_len(params$dimensionality), n.neighbors = params$umap_neighbors, min.dist = params$umap_min_dist, verbose = FALSE)

  if (return_markers) {
    # Compute markers only if requested, not during fitness evaluation.
    data_obj.markers <- Seurat::FindAllMarkers(seurat_obj, only.pos = TRUE)
    return(list(data_obj = seurat_obj, markers = data_obj.markers))
  }

  return(seurat_obj)
}