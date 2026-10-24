optimization_params <- c("scale_factor", "nfeatures", "dimensionality",
    "resolution", "umap_neighbors", "umap_min_dist")

# data loading - format check
load_count_data <- function(data_dir, gene_column) {
    h5_file <- list.files(data_dir, pattern = "\\.h5$", full.names = TRUE)[1]
    if (!is.na(h5_file)) {
        # check if hdf5r is installed to read in the file
        if (!requireNamespace("hdf5r", quietly = TRUE)) {
            stop(
                "The package 'hdf5r' is necessary to use the .h5 file type.
                Install it with install.packages('hdf5r').",
                call. = FALSE
            )
        }
        counts <- Seurat::Read10X_h5(filename = h5_file)
    } else {
        counts <- Seurat::Read10X(data.dir = data_dir, gene.column =
            gene_column)
    }

    # replace underscores with dashes in the feature names
    rownames(counts) <- gsub("_", "-", rownames(counts))
    return(counts)
}

# checks for a successful write access in the given directory
check_write_permissions <- function(directory) {
    test_file_path <- file.path(directory, "permission-test")
    writeLines("test", con = test_file_path)
    if (file.exists(test_file_path)) {
        file.remove(test_file_path)
    }
}

# checks whether all parameters have valid data types and values
check_parameter_validity <- function(gene_column, max_iterations, mu, lambda,
    survival_strategy, probability_mutation, adapt_mutation_rate,
    scale_factor_min, nfeatures_min, dimensionality_min, resolution_min,
    umap_neighbors_min, umap_min_dist_min, scale_factor_max, nfeatures_max,
    dimensionality_max, resolution_max, umap_neighbors_max, umap_min_dist_max,
    num_processes, qc_nfeature_min, qc_nfeature_max, qc_percent_mt_max,
    mt_gene_pattern, qc_gene_min_cells) {
    checkmate::assert_int(gene_column, lower = 1)
    checkmate::assert_int(max_iterations, lower = 1)
    checkmate::assert_int(mu, lower = 1)
    checkmate::assert_int(lambda, lower = 1)
    checkmate::assert_choice(survival_strategy, c("plus", "comma"))
    checkmate::assert_number(probability_mutation, lower = 0, upper = 1)
    checkmate::assert_number(adapt_mutation_rate, lower = 0, upper = 1)

    checkmate::assert_number(scale_factor_min, lower = 1000)
    checkmate::assert_number(scale_factor_max, lower = scale_factor_min,
        upper = 50000)

    checkmate::assert_number(nfeatures_min, lower = 500)
    checkmate::assert_number(nfeatures_max, lower = nfeatures_min, upper = 5000)

    checkmate::assert_number(dimensionality_min, lower = 5)
    checkmate::assert_number(dimensionality_max, lower = dimensionality_min,
        upper = 50)

    checkmate::assert_number(resolution_min, lower = 0.1)
    checkmate::assert_number(resolution_max, lower = resolution_min, upper = 2)

    checkmate::assert_number(umap_neighbors_min, lower = 5)
    checkmate::assert_number(umap_neighbors_max, lower = umap_neighbors_min,
        upper = 50)

    checkmate::assert_number(umap_min_dist_min, lower = 0.01)
    checkmate::assert_number(umap_min_dist_max, lower = umap_min_dist_min,
        upper = 0.7)

    checkmate::assert_int(num_processes, lower = 1,
        parallelly::availableCores())

    checkmate::assert_int(qc_nfeature_min, lower = 0)
    checkmate::assert_int(qc_nfeature_max, lower = qc_nfeature_min + 2)
    checkmate::assert_number(qc_percent_mt_max, lower = 0, upper = 100)
    checkmate::assert_string(mt_gene_pattern, min.chars = 1)
    checkmate::assert_int(qc_gene_min_cells, lower = 0)
}

# creates a sc_clustering configuration with the given values
create_sc_clustering_config <- function(scale_factor, nfeatures,
    dimensionality, resolution, umap_neighbors, umap_min_dist) {
    return(list(
        scale_factor = scale_factor,
        nfeatures = nfeatures,
        dimensionality = dimensionality,
        resolution = resolution,
        umap_neighbors = umap_neighbors,
        umap_min_dist = umap_min_dist
    ))
}

# creates a configuration as named vector with the given values
create_config <- function(results_dir, gene_column, max_iterations, mu,
    lambda, survival_strategy, probability_mutation, adapt_mutation_rate,
    scale_factor_min, nfeatures_min, dimensionality_min, resolution_min,
    umap_neighbors_min, umap_min_dist_min, scale_factor_max, nfeatures_max,
    dimensionality_max, resolution_max, umap_neighbors_max, umap_min_dist_max,
    qc_nfeature_min, qc_nfeature_max, qc_percent_mt_max, mt_gene_pattern,
    qc_gene_min_cells) {
    return(list(evolution = list(
            log = list(population = TRUE, stats = c("worst", "mean", "best")),
            gene_column = gene_column,
            init_default = TRUE,
            max_iterations = max_iterations,
            mu = mu,
            lambda = lambda,
            optimization_params = optimization_params,
            survival_strategy = survival_strategy,
            fitness_func = "Calinski_Harabasz",
            probability_mutation = probability_mutation,
            mutation_type = "adaptive",
            adapt_mutation_rate = adapt_mutation_rate),
        
        sc_clustering = list(
            default = create_sc_clustering_config(10000, 2000, 10, 0.5,
                30, 0.3),
            lower = create_sc_clustering_config(
                scale_factor_min, nfeatures_min, dimensionality_min,
                resolution_min, umap_neighbors_min, umap_min_dist_min),
            upper = create_sc_clustering_config(
                scale_factor_max, nfeatures_max, dimensionality_max,
                resolution_max, umap_neighbors_max, umap_min_dist_max),
            quality_thresholds = list(
                nFeature_RNA_min = qc_nfeature_min,
                nFeature_RNA_max = qc_nfeature_max,
                percent_mt = qc_percent_mt_max,
                mt_gene_pattern = mt_gene_pattern,
                gene_min_cells = qc_gene_min_cells)),
        
        results = list(
            dir_name = results_dir,
            file_scatter_fitness = "scatter-fitness.png",
            file_seurat_params = "seurat-params.json",
            file_umap_coordinates = "umap-coordinates.tsv",
            file_cluster_identities = "cluster-identities.tsv",
            file_umap_plot = "umap-static-default-seurat.png",
            file_umap_plot_interactive = "umap-interactive.html",
            dir_umap_plot_interactive = "umap-interactive_files",
            file_umap_plot_static = "umap-static.png")))
}

internal_run_optimization <- function(config, data_dir, results_dir) {
    count_data <- load_count_data(data_dir, config$evolution$gene_column)

    # exchange data with socket workers
    parallelMap::parallelExport(objnames = c("get_umap_coords",
        "fitness_function", "sc_clustering_pipeline", "config",
        "count_data", "check_is_named_list_or_matrix", "select_cells"))

    evolution_result <- run_evolution(config, count_data)
    best_candidate <- extract_best_candidate(count_data, evolution_result,
        config)
    message("Optimization finished successfully")

    message("Writing result files...")
    write_result(count_data = count_data, res = c(evolution_result,
        best_candidate), results_dir = results_dir, config = config,
        fitness_fun = fitness_function_factory(
            function_selection = config$evolution$fitness_func,
            count_data = count_data, config = config))

    names(best_candidate$best.params) <- optimization_params
    return(best_candidate$best.params)
}

#' Optimizes UMAP-based clustering pipeline parameters for single-cell
#' RNA-seq data using an evolutionary algorithm.
#'
#' This function optimizes UMAP-based clustering pipeline parameters for
#' single-cell RNA-seq data. Therefore, it treats the Seurat preprocessing and
#' clustering workflow as an optimization problem, identifying the most
#' appropriate configuration by the Calinski-Harabasz index.
#' It evaluates clustering quality by measuring the ratio of between-cluster
#' dispersion to within-cluster dispersion, given a set of data points in UMAP
#' space and their cluster labels, which are derived from the PCA space.
#'
#' The core idea is:
#'
#' 1. Load the 10X or h5 single-cell dataset.
#' 2. Run a fixed Seurat workflow with a small set of tunable parameters.
#' 3. Extract UMAP coordinates and cluster identities.
#' 4. Score the resulting clustering with the Calinski-Harabasz index.
#' 5. Use an evolutionary algorithm to explore the parameter space.
#' 6. Save the search result and export reporting artifacts for downstream
#' analysis.
#'
#' @section Parameter Optimization:
#'
#' The optimization includes the following parameters:
#'
#' - `scale_factor`
#' - `nfeatures`
#' - `dimensionality`
#' - `resolution`
#' - `umap_neighbors`
#' - `umap_min_dist`
#'
#' Each parameter is optimized within predefined upper and lower bounds which
#' can be set with the `[parameter]_min` and `[parameter]_max` function
#' parameters.
#' If `[parameter]_min = [parameter]_max`, the parameter does not evolve.
#' Note that the default parameters are always evaluated nonetheless.
#' The default parameter values are based on the [Seurat PBMC 3K tutorial](
#' https://satijalab.org/seurat/articles/pbmc3k_tutorial).
#' Although all six parameters are numerically encoded as floating-point
#' values during optimization, some are interpreted as integers later:
#' `dimensionality` is rounded before use and `nfeatures` and
#' `umap_neighbors` for analysis views.
#'
#' @section Evolution Configuration:
#'
#' The evolution and clustering behavior can be modified by the function
#' parameters:
#'
#' - which parameters are fixed and therefore excluded from exploration
#' - budget and population sizes for the evolutionary algorithm
#' - parent and survival selectors
#' - mutation regime

#' The speed of mutation-rate adaptation is controlled by the
#' `adapt_mutation_rate` parameter.
#'
#' The **parent and survival selection** uses a tournament selection for
#' mating and a greedy replacement for survival.
#' The survival strategy itself is configurable with the `survival_strategy`
#' parameter:
#'
#' - `plus` means `(mu + lambda)` survival
#' - `comma` means `(mu, lambda)` survival
#'
#' The default is `plus`, meaning parents and offspring compete together for
#' survival.
#'
#' @section Dataset:
#'
#' The dataset to run the optimization with may be an h5 file or a 10X
#' Genomics dataset.
#' Both types are accepted by Seurat and read in by `Read10X(...)` /
#' `Read10X_h5(...)` and `CreateSeuratObject(...)`.
#' After loading of the dataset, quality control filtering is applied. Its
#' thresholds and mitochondrial gene-name pattern can be set with the `qc_*`
#' and `mt_gene_pattern` function parameters:
#'
#' - mitochondrial percentage is computed using `mt_gene_pattern`
#' - genes detected in fewer than `qc_gene_min_cells` cells are excluded
#' - cells are retained only if:
#' - `nFeature_RNA > qc_nfeature_min`
#' - `nFeature_RNA < qc_nfeature_max`
#' - `percent_mt < qc_percent_mt_max`
#'
#' The defaults are `qc_nfeature_min = 200`, `qc_nfeature_max = 2500`,
#' `qc_percent_mt_max = 5`, `mt_gene_pattern = "^MT-"`, and
#' `qc_gene_min_cells = 3`.
#'
#' @section Seurat Pipeline Execution:
#'
#' For a candidate parameter set, the routine performs the following sequence:
#'
#' 1. Normalize counts using the candidate `scale_factor`
#' 2. Select variable features using candidate `nfeatures`
#' 3. Scale all genes
#' 4. Run PCA
#' 5. Build the neighbor graph using the first `dimensionality` PCs
#' 6. Cluster cells using candidate `resolution`
#' 7. Compute UMAP using candidate `umap_neighbors` and `umap_min_dist`
#'
#' This is an important design decision: the search is not optimizing a
#' custom embedding model directly. It is optimizing
#' hyperparameters of a standard Seurat workflow, and all fitness is derived
#' from the resulting UMAP-plus-clustering state.
#'
#' @section Parallelization:
#'
#' The evolution uses `parallelMap`, at level `ecr.evaluateFitness` to
#' **parallelize** the fitness evaluation for all individuals of the population.
#' Therefore, each individual can be evaluated in a different parallelization
#' unit. This is recommended because evaluating one individual means running a
#' full Seurat pipeline, which is computationally expensive.
#' The number of processes can be set with the `num_processes` function
#' parameter.
#' You can then disable parallelization by setting `num_processes = 1`.
#'
#' @param data_dir the directory with a .h5 file or 10X Genomics dataset
#' (expected files: barcodes.tsv, genes.tsv, matrix.mtx) to run the
#' optimization with
#' @param results_dir the directory to save the result files of the
#' optimization, the directory is created if it does not exist
#' @param gene_column the position of the gene names in the genes.tsv file on
#' a 10x Genomics dataset
#' @param max_iterations the maximum number of iterations of the evolution
#' @param mu (greek letter \eqn{\mu}) the number of individuals in the
#' population (commonly referred to as the parent population size)
#' @param lambda the number of individuals generated in each generation
#' @param survival_strategy the survivor-selection strategy: "plus" (mu +
#' lambda) selects the next population of `mu` individuals from the parents and
#' offspring combined, so parents may survive; "comma" (mu, lambda) selects the
#' next population only from the `lambda` offspring, so parents do not survive
#' directly. For "comma", `lambda` must be at least `mu`.
#' @param probability_mutation the probability that a mutation operator will
#' be applied to a child - refers only to the application of the mutation
#' operator, not to the probability of mutating individual genes of the
#' respective child
#' @param adapt_mutation_rate controls how fast mutation rates can evolve
#' @param scale_factor_min the lower bound for the scale_factor parameter in
#' the optimization
#' @param nfeatures_min the lower bound for the nfeatures parameter in the
#' optimization
#' @param dimensionality_min the lower bound for the dimensionality parameter
#' in the optimization
#' @param resolution_min the lower bound for the resolution parameter in
#' the optimization
#' @param umap_neighbors_min the lower bound for the umap_neighbors parameter
#' in the optimization
#' @param umap_min_dist_min the lower bound for the umap_min_dist parameter
#' in the optimization
#' @param scale_factor_max the upper bound for the scale_factor parameter in
#' the optimization
#' @param nfeatures_max the upper bound for the nfeatures parameter in the
#' optimization
#' @param dimensionality_max the upper bound for the dimensionality parameter
#' in the optimization
#' @param resolution_max the upper bound for the resolution parameter in
#' the optimization
#' @param umap_neighbors_max the upper bound for the umap_neighbors parameter
#' in the optimization
#' @param umap_min_dist_max the upper bound for the umap_min_dist parameter
#' in the optimization
#' @param num_processes the number of processes to run the fitness evaluation
#' with - default is the minimum of number of individuals and available cores
#' @param qc_nfeature_min minimum detected features for retaining a cell
#' @param qc_nfeature_max maximum detected features for retaining a cell
#' @param qc_percent_mt_max maximum mitochondrial percentage for retaining a
#' cell
#' @param mt_gene_pattern regular expression identifying mitochondrial genes
#' @param qc_gene_min_cells minimum number of cells in which a gene must be
#' detected to be retained
#' @export
#' @returns the optimized parameters of the evolution as named vector
#' @examples
#' # load package example data
#' data_dir <- system.file("extdata", "exampledata", package="geneticumap")
#'
#' # run simple example with multiple processes
#' results_dir_1 <- file.path(tempdir(), "geneticumap-example-1")
#' run_optimization(data_dir = data_dir, results_dir = results_dir_1, mu = 2,
#'     lambda = 4, max_iterations = 8)
#'
#' # run simple example with single process and fix dimensionality
#' results_dir_2 <- file.path(tempdir(), "geneticumap-example-2")
#' run_optimization(data_dir = data_dir, results_dir = results_dir_2,
#'     dimensionality_min = 6, dimensionality_max = 6, num_processes = 1)
run_optimization <- function(data_dir, results_dir = file.path(data_dir,
    "results"), gene_column = 1, max_iterations = 10, mu = 5, lambda = 10,
    survival_strategy = "plus", probability_mutation = 1,
    adapt_mutation_rate = 0.7, scale_factor_min = 1000, nfeatures_min = 500,
    dimensionality_min = 5, resolution_min = 0.1, umap_neighbors_min = 5,
    umap_min_dist_min = 0.01, scale_factor_max = 50000, nfeatures_max = 5000,
    dimensionality_max = 50, resolution_max = 2, umap_neighbors_max = 50,
    umap_min_dist_max = 0.7,
    num_processes = min(lambda, parallelly::availableCores()),
    qc_nfeature_min = 200, qc_nfeature_max = 2500,
    qc_percent_mt_max = 5, mt_gene_pattern = "^MT-",
    qc_gene_min_cells = 3) {
    check_parameter_validity(gene_column, max_iterations, mu, lambda,
        survival_strategy, probability_mutation, adapt_mutation_rate,
        scale_factor_min, nfeatures_min, dimensionality_min, resolution_min,
        umap_neighbors_min, umap_min_dist_min, scale_factor_max, nfeatures_max,
        dimensionality_max, resolution_max, umap_neighbors_max,
        umap_min_dist_max, num_processes, qc_nfeature_min, qc_nfeature_max,
        qc_percent_mt_max, mt_gene_pattern, qc_gene_min_cells)

    message("Initializing optimization...")

    config <- create_config(results_dir, gene_column, max_iterations, mu,
        lambda, survival_strategy, probability_mutation, adapt_mutation_rate,
        scale_factor_min, nfeatures_min, dimensionality_min, resolution_min,
        umap_neighbors_min, umap_min_dist_min, scale_factor_max, nfeatures_max,
        dimensionality_max, resolution_max, umap_neighbors_max,
        umap_min_dist_max, qc_nfeature_min, qc_nfeature_max,
        qc_percent_mt_max, mt_gene_pattern, qc_gene_min_cells)

    # create results directory
    if (!dir.exists(results_dir)) dir.create(results_dir, showWarnings = TRUE)
    check_write_permissions(results_dir)

    # setup parallelization with socket workers
    on.exit(expr = parallelMap::parallelStop(), add = TRUE)
    parallelMap::parallelRegisterLevels(package = "ecr", levels =
        "evaluateFitness")
    parallelMap::parallelLibrary("clusterCrit", "ecr", "BBmisc", "Seurat",
        level = "ecr.evaluateFitness")
    parallelMap::parallelStart(mode = "socket", cpus = num_processes)

    result_params <- internal_run_optimization(config, data_dir, results_dir)

    # stop socket workers
    parallelMap::parallelStop()

    return(result_params)
}