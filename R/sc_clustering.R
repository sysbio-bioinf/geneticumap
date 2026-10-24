options(Seurat.warn.umap.uwot = FALSE)

check_is_named_list_or_matrix <- function(count_data) {
    checkmate::assert(checkmate::check_list(count_data, names = "unique"),
        checkmate::check_class(count_data, "dgCMatrix"),
        combine = "or", .var.name = "count_data")
}

select_cells <- function(seurat_obj, config) {
    seurat_meta <- seurat_obj[[]]
    keep_cells <- seurat_meta$nFeature_RNA > config$sc_clustering$
        quality_thresholds$nFeature_RNA_min & seurat_meta$nFeature_RNA <
        config$sc_clustering$quality_thresholds$nFeature_RNA_max &
        seurat_meta$percent_mt < config$sc_clustering$quality_thresholds$
        percent_mt
    seurat_obj <- seurat_obj[, keep_cells]
    return(seurat_obj)
}

sc_clustering_pipeline <- function(count_data, config, params, return_markers
    = FALSE) {
    check_is_named_list_or_matrix(count_data)

    seurat_obj <- Seurat::CreateSeuratObject(counts = count_data, project =
        "geneticumap", min.cells = config$sc_clustering$quality_thresholds$
        gene_min_cells, min.features = config$sc_clustering$
        quality_thresholds$nFeature_RNA_min)

    # quality control and cell selection
    seurat_obj[["percent_mt"]] <- Seurat::PercentageFeatureSet(seurat_obj,
        pattern = config$sc_clustering$quality_thresholds$mt_gene_pattern)
    seurat_obj <- select_cells(seurat_obj, config)

    # round integer-valued parameters (optimizer works with continuous floats)
    params$nfeatures <- round(params$nfeatures)
    params$dimensionality <- round(params$dimensionality)
    params$umap_neighbors <- round(params$umap_neighbors)

    seurat_obj <- Seurat::NormalizeData(seurat_obj, normalization.method =
        "LogNormalize", scale.factor = params$scale_factor, verbose = FALSE)

    # identification of highly variable features
    seurat_obj <- Seurat::FindVariableFeatures(seurat_obj, selection.method =
        "vst", nfeatures = params$nfeatures, verbose = FALSE)

    # linear transformation (scaling)
    all.genes <- rownames(seurat_obj)
    seurat_obj <- Seurat::ScaleData(seurat_obj, features = all.genes, verbose
        = FALSE)

    seurat_obj <- Seurat::RunPCA(seurat_obj, features =
        Seurat::VariableFeatures(object = seurat_obj), verbose = FALSE)

    # cluster the cells
    seurat_obj <- Seurat::FindNeighbors(seurat_obj, dims =
        seq_len(params$dimensionality), verbose = FALSE)
    seurat_obj <- Seurat::FindClusters(seurat_obj, resolution =
        params$resolution, verbose = FALSE)

    seurat_obj <- Seurat::RunUMAP(seurat_obj, dims =
        seq_len(params$dimensionality), n.neighbors = params$umap_neighbors,
        min.dist = params$umap_min_dist, verbose = FALSE)

    if (return_markers) { # Compute markers if requested
        data_obj.markers <- Seurat::FindAllMarkers(seurat_obj, only.pos = TRUE)
        return(list(data_obj = seurat_obj, markers = data_obj.markers))
    }
    return(seurat_obj)
}
