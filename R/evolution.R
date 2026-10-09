get_umap_coords <- function(count_data, config, parm) {
  params <- as.list(parm)
  names(params) <- config$evolution$optimization_params
  seurat_obj <- sc_clustering_pipeline(count_data = count_data, config = config, params = params)
  SeuratObject::FetchData(seurat_obj, c("umap_1", "umap_2", "ident"))
}

fitness_function_factory <- function(function_selection, count_data, config) {
  return(\(umap_coors) fitness_function(function_selection = function_selection, umap_coordinates = umap_coors, count_data = count_data, config = config))
}

#fitness function definition
fitness_function <- function(function_selection, umap_coordinates, count_data, config) {
  umap_coors <- tryCatch(
    get_umap_coords(count_data = count_data, config = config, parm = umap_coordinates),
    error = function(e) {
      warning(sprintf("Fitness evaluation failed for candidate: %s", e$message), call. = FALSE)
      return(NULL)
    }
  )

  if (is.null(umap_coors)) {
    if (function_selection == "Calinski_Harabasz") {
      return(-1e9)
    }
    return(Inf)
  }

  umap_matrix <- cbind(umap_coors$umap_1, umap_coors$umap_2)
  unlist(clusterCrit::intCriteria(umap_matrix, as.integer(umap_coors$ident), function_selection))
}

fitness_function_is_minimize <- function(function_selection) {
  if (function_selection == "Calinski_Harabasz")
    return(FALSE)
  else
    return(TRUE)
}

#helper functions for mutation
#' @exportS3Method NULL
mutate.gauss <- function(X, p = 1, sdev = 0.01, lower, upper) {
  n <- length(X)
  mut.idx <- stats::runif(n) < p
  sdev[is.na(sdev)] <- 0 #set all NA values to 0
  mut.noise <- stats::rnorm(sum(mut.idx), mean = 0, sd = sdev)
  X[mut.idx] <- X[mut.idx] + mut.noise
  pmin(pmax(lower, X), upper)
}

check_numeric_mutator_arguments <- function(ind, lower, upper) {
  checkmate::assertNumeric(x = ind, min.len = 1L, any.missing = FALSE, all.missing = FALSE)
  checkmate::assertNumeric(x = lower, any.missing = FALSE, all.missing = FALSE)
  checkmate::assertNumeric(x = upper, any.missing = FALSE, all.missing = FALSE)
}

#mutation operator setup
setup_mutator <- function(mutation_type, adapt_mutation_rate, n_dim) {
  return(if (mutation_type == "fixedGauss") {
    ecr::makeMutator(
      function(ind, p = 1L, sdev = 0.05, lower, upper, ...) {
        checkmate::assertNumber(p, lower = 0, finite = TRUE, na.ok = FALSE)
        checkmate::assertNumber(sdev, lower = 0, finite = TRUE, na.ok = FALSE)

        mutate.gauss(ind, p, sdev, lower, upper)
      },
      supported = "float"
    )
  } else if (mutation_type == "annealing") {
    ecr::makeMutator(
      function(ind, p = 1L, sdev = 0.05, n.iter = 1, lower, upper, ...) {
        checkmate::assertNumber(p, lower = 0, finite = TRUE, na.ok = FALSE)

        mutate.gauss(ind, p, sdev, lower, upper)
      },
      supported = "float"
    )
  } else if (mutation_type == "adaptive") {
    ecr::makeMutator(
      function(ind, p = 1L, sdev = 0.05, lower, upper, ...) {
        check_numeric_mutator_arguments(ind, lower, upper)
        checkmate::assertNumber(p, lower = 0, finite = TRUE, na.ok = FALSE)
        checkmate::assertNumber(sdev, lower = 0, finite = TRUE, na.ok = FALSE)
        if (length(ind) %% 2 != 0) stop("Number of mutation rates does not match number of parameters!")

        n.parm <- n_dim / 2
        params <- ind[1:n.parm]
        mutation.rates <- ind[(n.parm + 1):length(ind)]

        mutation.rates <- mutation.rates * (stats::rlnorm(n.parm)^adapt_mutation_rate)

        range.size <- upper[1:n.parm] - lower[1:n.parm]
        sdev <- range.size * mutation.rates
        params <- mutate.gauss(params, p, sdev, lower[1:n.parm], upper[1:n.parm])

        return(c(params, mutation.rates))
      },
      supported = "float"
    )
  })
}

run_evolution <- function(config, count_data) {
  n_dim <- if (config$evolution$mutation_type == "adaptive") 2L * length(config$evolution$optimization_params) else length(config$evolution$optimization_params)

  if (config$evolution$init_default) {
    default_params <- unlist(config$sc_clustering$default, use.names = FALSE)
    if (config$evolution$mutation_type == "adaptive") {
      default_mut_rates <- rep(0.1, length(config$evolution$optimization_params))
      initial_solutions <- list(c(default_params, default_mut_rates))
    } else {
      initial_solutions <- list(default_params)
    }
  } else {
    initial_solutions <- NULL
  }

  lower <- unlist(sapply(config$evolution$optimization_params, \(parm) config$sc_clustering$lower[[parm]]), use.names = FALSE)
  upper <- unlist(sapply(config$evolution$optimization_params, \(parm) config$sc_clustering$upper[[parm]]), use.names = FALSE)
  #adaptive bounds for mutation rates if using adaptive mutation
  if (config$evolution$mutation_type == "adaptive") {
    lower <- c(lower, rep.int(0.00005, length(config$evolution$optimization_params)))
    upper <- c(upper, rep.int(1, length(config$evolution$optimization_params)))
  }

  #run single-objective evolutionary algorithm
  result <- evolution_loop(
    fitness.fun = fitness_function_factory(function_selection = config$evolution$fitness_func, config = config, count_data = count_data),
    fitness.fun.name = config$evolution$fitness_func,
    minimize = fitness_function_is_minimize(config$evolution$fitness_func),
    n.objectives = 1L,
    representation = "float",
    n.dim = n_dim,
    lower = lower,
    upper = upper,
    mu = config$evolution$mu,
    lambda = config$evolution$lambda,
    p.recomb = config$evolution$probability_recombination,
    p.mut = config$evolution$probability_mutation,
    survival.strategy = config$evolution$survival_strategy,
    initial.solutions = initial_solutions,
    max.iter = config$evolution$max_iterations,
    mutator = setup_mutator(mutation_type = config$evolution$mutation_type, adapt_mutation_rate = config$evolution$adapt_mutation_rate, n_dim = n_dim),
    parent.selector = ecr::selTournament, #parent selector for single-objective task
    survival.selector = ecr::selGreedy, #survival selector for single-objective task
    log.stats = config$evolution$log$stats,
    log.pop = config$evolution$log$population
  )
  return(result)
}

write_result <- function(count_data, res, results_dir, fitness_fun, config) {
  # add marker genes, umap coordinates and cluster identities of the best candidate to results object
  # extract and round parameters of best candidate
  best.params <- as.list(res$best.x[[1]][seq_along(config$evolution$optimization_params)])
  names(best.params) <- config$evolution$optimization_params
  digits <- c(scale_factor = 2L, nfeatures = 0L, dimensionality = 0L, resolution = 2L, umap_neighbors = 0L, umap_min_dist = 2L)
  best.params[names(digits)] <- Map(\(x, d) if (d == 0L) as.integer(round(x)) else round(as.numeric(x), d), best.params[names(digits)], digits)
  # compute umap and find markers for best parameters
  best.result <- sc_clustering_pipeline(count_data = count_data, config = config, params = best.params, return_markers = TRUE)
  res$best.markers <- best.result$markers
  res$best.umap <- Seurat::Embeddings(best.result$data_obj, "umap")
  res$best.clusters <- best.result$data_obj$seurat_clusters
  res$best.seurat <- best.result$data_obj

  # export parameters, umap coordinates and cluster identities of best candidate
  seurat_params_file <- file.path(results_dir, config$results$file_seurat_params)
  jsonlite::write_json(best.params, seurat_params_file, pretty = TRUE, auto_unbox = TRUE)

  umap_coordinates_file <- file.path(results_dir, config$results$file_umap_coordinates)
  utils::write.table(data.frame(cell_id = rownames(res$best.umap), res$best.umap, row.names = NULL, check.names = FALSE), umap_coordinates_file, sep = "\t", quote = FALSE, row.names = FALSE)

  cluster_identities_file <- file.path(results_dir, config$results$file_cluster_identities)
  utils::write.table(data.frame(cell_id = names(res$best.clusters), cluster = as.character(res$best.clusters), row.names = NULL, check.names = FALSE), cluster_identities_file, sep = "\t", quote = FALSE, row.names = FALSE)

  # save scatter_fitness_plot and interactive and static umap plot of best candidate + umap with default seurat params for reference
  scatter_fitness_file <- file.path(results_dir, config$results$file_scatter_fitness)
  scatter_plot_from_res(fitness_fun = fitness_fun, config = config, output_file = scatter_fitness_file, res = res)

  umap_plot_file <- file.path(results_dir, config$results$file_umap_plot)
  plot_umap(count_data = count_data, config = config, ind = unlist(config$sc_clustering$default, use.names = FALSE), filename = umap_plot_file)

  umap_plot_html <- file.path(results_dir, config$results$file_umap_plot_interactive)
  umap_plot_static <- file.path(results_dir, config$results$file_umap_plot_static)

  # check for optional package htmlwidgets
  if (requireNamespace("htmlwidgets", quietly = TRUE)) {
    plot_interactive_umap_from_res(output_file_html = umap_plot_html, output_file_static = umap_plot_static, res = res)
  } else {
    message("Optional package 'htmlwidgets' is not installed. No interactive umap plot is generated.")
  }

  return(best.params)
}