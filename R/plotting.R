library(ggplot2, quietly = TRUE)

plot_umap <- function(count_data, config, ind, filename, ...) {
    ind <- as.list(ind)
    names(ind) <- config$evolution$optimization_params
    umap_plot <- Seurat::DimPlot(sc_clustering_pipeline(count_data = count_data, config = config, params = ind), reduction = "umap", ...) +
        ggplot2::scale_color_discrete(name = "Cluster", labels = function(x) paste("Cluster", x)) +
        ggplot2::theme(plot.margin = ggplot2::margin(5, 5, 0, 5))
    ggplot2::ggsave(filename = filename, plot = umap_plot, width = 210, height = 170, units = "mm", dpi = 600)
}

plot_interactive_umap_from_res <- function(output_file_html, output_file_static, res) {
    # extract coordinates, cluster idents and markers from results object
    umap <- data.frame(cell_id = rownames(res$best.umap), res$best.umap, cluster = as.character(res$best.clusters[rownames(res$best.umap)]), check.names = FALSE)
    # ensure legend orders clusters numerically rather than alphabetically
    cluster_ids <- sort(unique(as.integer(umap$cluster)))
    cluster_labels <- paste("Cluster", cluster_ids)
    cluster_colors <- stats::setNames(grDevices::hcl.colors(length(cluster_ids), palette = "Dynamic"), cluster_labels)
    umap$cluster_label <- factor(paste("Cluster", umap$cluster), levels = cluster_labels)
    markers <- res$best.markers; if (is.null(markers$gene)) markers$gene <- rownames(markers)

    # design hover text
    top20 <- lapply(split(markers, as.character(markers$cluster)), \(df) utils::head(df$gene[order(df$avg_log2FC, decreasing = TRUE)], 20))
    hover <- stats::setNames(vapply(names(top20), \(cl) {
        genes <- top20[[cl]]
        gene_text <- paste(sapply(seq(1, length(genes), 5), \(i) paste(genes[i:min(i + 4, length(genes))], collapse = ", ")), collapse = "<br>")
        paste0("Cluster: ", cl, "<br>Top 20 markers:<br>", gene_text)
    }, character(1)), names(top20))

    # generate and save plot
    p <- plotly::plot_ly(umap, x = ~umap_1, y = ~umap_2, type = "scattergl", mode = "markers", color = ~cluster_label, colors = unname(cluster_colors), text = ~hover[cluster], hoverinfo = "text")
    htmlwidgets::saveWidget(p, output_file_html, selfcontained = TRUE, libdir = NULL)
    static_colors <- stats::setNames(unname(cluster_colors), as.character(cluster_ids))
    umap_plot <- Seurat::DimPlot(res$best.seurat, reduction = "umap") +
        ggplot2::scale_color_manual(name = "Cluster", values = static_colors, breaks = names(static_colors), labels = cluster_labels, drop = FALSE) +
        ggplot2::theme(plot.margin = ggplot2::margin(5, 5, 0, 5))
    ggplot2::ggsave(filename = output_file_static, plot = umap_plot, width = 210, height = 170, units = "mm", dpi = 600)
    invisible(p)
}

# x axis = generation, y axis = fitness of each individual per generation, blue dashed line = fitness of default parameters
scatter_singledim <- function(ecr_result, file_name, ymax, default_chi) {
    population <- ecr_result$log$env$pop
    indexed_fitness <- as.data.frame(do.call(rbind, lapply(seq_along(population), function(i) {
        t(vapply(population[[i]]$fitness, function(x) c(x, i - 1), FUN.VALUE = numeric(2)))
    })))
    colnames(indexed_fitness) <- c("Fitness", "Generation")
    scatter_plot <- ggplot2::ggplot(indexed_fitness, ggplot2::aes(x = .data[["Generation"]], y = .data[["Fitness"]])) +
        ggplot2::coord_cartesian(ylim = c(0, ymax)) +
        ggplot2::geom_point(size = 0.3)
    ggplot2::ggsave(filename = file_name, plot = scatter_plot, width = 210, height = 170, units = "mm", dpi = 600)
}

calc_default_chi <- function(fitness_fun, config) {
    p <- unlist(config$sc_clustering$default, use.names = FALSE)
    as.numeric(fitness_fun(p))
}

calc_ymax <- function(ecr_result, default_chi, pad = 0.05) {
    pop <- ecr_result$log$env$pop
    vals <- unlist(lapply(pop, \(g) as.numeric(g$fitness)))
    max(c(vals, default_chi), na.rm = TRUE) * (1 + pad)
}

# wrapper with dynamic ymax and default_chi based on results object
scatter_plot_from_res <- function(fitness_fun, config, output_file, res) {
    default_chi <- calc_default_chi(fitness_fun = fitness_fun, config = config)
    ymax <- calc_ymax(res, default_chi)
    scatter_singledim(res, output_file, ymax, default_chi)
}
