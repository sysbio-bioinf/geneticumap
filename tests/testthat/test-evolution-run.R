# delete existing test directory, return directory path to write results to
prepare_test_directory <- function() {
  # use temporary results directory
  results_dir <- file.path(tempdir(), "test")
  message("Temporary diectory for test results: ", results_dir)

  # delete directory to remove al files
  unlink(results_dir, recursive = TRUE)

  return(results_dir)
}

# check all result files are written to results directory and
# seurat-parameters.json contains all optimized parameters
check_result_files <- function(results_dir) {
  expected_files <- list(
      file_scatter_fitness = "scatter-fitness.png",
      file_seurat_params = "seurat-params.json",
      file_umap_coordinates = "umap-coordinates.tsv",
      file_cluster_identities = "cluster-identities.tsv",
      file_umap_plot = "umap-static-default-seurat.png",
      file_umap_plot_interactive = "umap-interactive.html",
      file_umap_plot_static = "umap-static.png"
    )
  result_files <- list.files(results_dir)

  expect_true(all(expected_files %in% result_files))

  # test seurat-parameters.json is valid
  optimization_params <- c("scale_factor", "nfeatures", "dimensionality", "resolution", "umap_neighbors", "umap_min_dist")
  file_path_seurat_params <- file.path(results_dir, expected_files$file_seurat_params)
  library(jsonlite)
  seurat_params <- jsonlite::fromJSON(file_path_seurat_params, simplifyVector = FALSE)
  missing_params <- setdiff(optimization_params, names(seurat_params))
  expect_equal(length(missing_params), 0)
}

check_scale_factor_in_range <- function(scale_factor, min = 1000, max = 50000) {
  expect_gte(scale_factor, min)
  expect_lte(scale_factor, max)
}

check_nfeatures_in_range <- function(scale_factor, min = 500, max = 5000) {
  expect_gte(scale_factor, min)
  expect_lte(scale_factor, max)
}

check_dimensionality_in_range <- function(scale_factor, min = 5, max = 50) {
  expect_gte(scale_factor, min)
  expect_lte(scale_factor, max)
}

check_resolution_in_range <- function(scale_factor, min = 0.1, max = 2) {
  expect_gte(scale_factor, min)
  expect_lte(scale_factor, max)
}

check_umap_neighbors_in_range <- function(scale_factor, min = 5, max = 50) {
  expect_gte(scale_factor, min)
  expect_lte(scale_factor, max)
}

check_umap_min_dist_in_range <- function(scale_factor, min = 0.01, max = 0.7) {
  expect_gte(scale_factor, min)
  expect_lte(scale_factor, max)
}

# run full evolution with default parameters
test_that("Evolution runs successfully (single-threaded)", {
  results_dir <- prepare_test_directory()

  # run evolution
  tryCatch({
    result_params <- run_optimization(data_dir = "../testdata", results_dir = results_dir, mu = 2, lambda = 4, max_iterations = 2, num_processes = 1)

    # check optimized parameters are in default ranges
    check_scale_factor_in_range(result_params$scale_factor)
    check_nfeatures_in_range(result_params$nfeatures)
    check_dimensionality_in_range(result_params$dimensionality)
    check_resolution_in_range(result_params$resolution)
    check_umap_neighbors_in_range(result_params$umap_neighbors)
    check_umap_min_dist_in_range(result_params$umap_min_dist)

    # check all expected result files are present
    check_result_files(results_dir)

  }, error = function(info) {
    fail(message = paste0("Test failed with error: ", info))
  })
})

# set invalid parameters and expect an error message
test_that("Parameter checks fail as expected", {
  # survival_strategy not in possible choices
  expect_error({
    run_optimization(data_dir = "../testdata", results_dir = results_dir, survival_strategy = "Something")
  })
  # dimensionality out of range (<50)
  expect_error({
    run_optimization(data_dir = "../testdata", results_dir = results_dir, survival_strategy = "Something")
  })
})

# set one optimized parameter to a fix value by setting _min = _max and expect the resulting parameter to stick with the fix value
test_that("Evolution runs with fix parameter (with socket workers)", {
  results_dir <- prepare_test_directory()

  # run evolution with fixed dimensionality
  fix_dimensionality <- 6
  default_dimensionality <- 10
  tryCatch({
    result_params <- run_optimization(data_dir = "../testdata", results_dir = results_dir, dimensionality_min = fix_dimensionality,
                                      dimensionality_max = fix_dimensionality, mu = 2, lambda = 4, max_iterations = 2)

    # check optimized parameters are in expected ranges
    check_scale_factor_in_range(result_params$scale_factor)
    check_nfeatures_in_range(result_params$nfeatures)
    expect_true(result_params$dimensionality == fix_dimensionality || result_params$dimensionality == default_dimensionality)
    check_resolution_in_range(result_params$resolution)
    check_umap_neighbors_in_range(result_params$umap_neighbors)
    check_umap_min_dist_in_range(result_params$umap_min_dist)

    # check all expected result files are present
    check_result_files(results_dir)

  }, error = function(info) {
    fail(message = paste0("Test failed with error: ", info))
  })
})
