# geneticUMAP

---

geneticUMAP optimizes UMAP-based clustering pipeline parameters for single-cell RNA-seq data, such as resolution, the number of principal components, and the number of neighbors. It uses an evolutionary algorithm to search for pipeline parameters for a given dataset. The Seurat preprocessing and clustering workflow is treated as an optimization problem, and the search identifies the configuration with the highest Calinski-Harabasz index.

The core idea is:

1. Load the 10X or h5 single-cell dataset.
2. Run a fixed Seurat workflow with a small set of tunable parameters.
3. Extract UMAP coordinates and cluster identities.
4. Score the resulting clustering with the Calinski-Harabasz index.
5. Use an evolutionary algorithm to explore the parameter space.
6. Save the search result and export reporting artifacts for downstream analysis.

The geneticUMAP package is developed at the [Institute of Medical Systems Biology of Ulm University](https://sysbio.uni-ulm.de).

---

## Associated Paper (Unpublished)

**geneticUMAP: Objective and Robust Parameter Optimization for scRNA-Seq Data Visualization**

**Authors:** Lisa-Maria Winter<sup>1,2</sup>, Hanna Wagner<sup>1</sup>, Dominik Dreiheller, Julian Schwab<sup>2</sup>, Johann M. Kraus<sup>1</sup>, and Hans A. Kestler<sup>1</sup>

<sup>1</sup> Ulm University, Institute of Medical Systems Biology  
<sup>2</sup> Boehringer Ingelheim Pharma GmbH & Co. KG

## Citation

The associated paper has not yet been published. Once published, please cite it, if it was helpful to you. Thank you.

## Web Application

For non-sensitive data you can use our [web application](https://fermat.informatik.uni-ulm.de/geneticumap) 
to run a parameter optimization in a comfortable user interface without any installation.

---

## Features

- **Automated Parameter Optimization:** Systematically tunes normalization scale factors, number of variable features, PCA dimensions, clustering resolution, UMAP neighbors, and minimum distance.
- **Self-Adaptive Evolutionary Search:** Utilizes adaptive mutation rates that evolve alongside parameters to navigate complex parameter landscapes efficiently.
- **Parallel Processing:** Integrated with `parallelMap` to evaluate candidate parameter configurations concurrently.
- **Comprehensive Artifact Generation:** Exports winner parameters as JSON, cluster and UMAP coordinate tables, static diagnostic plots, and interactive HTML UMAP visualizations with top cluster marker genes.

---

## Prerequisites & Installation

Make sure [R](https://www.r-project.org/) ist installed on your system.

The program has several dependencies on packages that need to be installed first.

You can install the package via Bioconductor:

```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("geneticumap")
```
---

## Quick Start

The primary entry point is the `run_optimization()` function. Below is a minimal example showing how to run an optimization on a local dataset:

```
library(geneticumap)

# Run optimization on a 10X or HDF5 dataset
result_parameters <- run_optimization(
  data_dir = "path/to/data",
  results_dir = "path/to/data/geneticumap_results",
  num_processes = 4
)

# Inspect the best identified parameter set
print(result_parameters)
```

## Usage

The `run_optimization()` function executes the optimization.
The function parameters let you configure:

- the folder paths for data and results
- the parameter bounds
- the iteration budget and population sizes
- the survival strategy and recombination and mutation probabilities
- the rate of mutation-rate adaptation
- the number of parallel processes

The fitness function is fixed to the Calinski-Harabasz index.
The default parameter set is included in the initial population, and the mutation regime is adaptive.

## Input Data

geneticUMAP accepts single-cell RNA-seq expression datasets supported by standard Seurat readers:

- 10X Genomics directory: A folder containing `matrix.mtx`, `barcodes.tsv`, and `features.tsv` (or `genes.tsv`) files. Loaded via `Seurat::Read10X()`.
- HDF5 file (.h5): A directory containing a standard Cell Ranger HDF5 output file. Loaded via `Seurat::Read10X_h5()`.

After loading, the pipeline applies dataset-agnostic quality control filtering:

- mitochondrial percentage is computed with pattern `^MT-`
- cells are retained only if:
    - `nFeature_RNA > 200`
    - `nFeature_RNA < 2500`
    - `percent_mt < 5`

These QC filters are hard-coded and apply to all datasets.

## Explored Parameters

The search operates on six tunable parameters:

- `scale_factor`
- `nfeatures`
- `dimensionality`
- `resolution`
- `umap_neighbors`
- `umap_min_dist`

The bounds for these parameters can be set with the corresponding function parameters to `run_optimization()`.

All six parameters are numerically encoded as floating-point values during optimization, but `dimensionality`, `nfeatures`, and `umap_neighbors` are rounded to integers before being passed to Seurat.

## Generated Results

The `run_optimization()` function produces the following result files in `results_dir`:

- `seurat-params.json`: Optimal pipeline parameter configuration identified by the search.
- `umap-coordinates.tsv`: Final UMAP coordinates (`umap_1`, `umap_2`) for all cells.
- `cluster-identities.tsv`: Cell barcode to cluster identity mappings.
- `scatter-fitness.png`: Diagnostic plot showing best, mean, and worst fitness across generations.
- `umap-static.png`: Static plot (`Seurat::DimPlot`) of the optimal UMAP embedding.
- `umap-static-default-seurat.png`: Static plot of the embedding produced with the default parameter set.
- `umap-interactive.html`: Self-contained interactive Plotly UMAP with hover tooltips displaying the top 20 cluster marker genes.

---

## Algorithmic Details

### Seurat Pipeline Execution

For each candidate parameter set, geneticUMAP executes a full standard Seurat pipeline:

1. Normalize counts with `Seurat::NormalizeData()` using the candidate `scale_factor`
2. Identify variable features with `Seurat::FindVariableFeatures()` using candidate `nfeatures`
3. Scale expression with `Seurat::ScaleData()`
4. Perform PCA with `Seurat::RunPCA()`
5. Compute neighbor graph with `Seurat::FindNeighbors()` using the first `dimensionality` PCs
6. Cluster cells with `Seurat::FindClusters()` using candidate `resolution`
7. Generate UMAP embedding with `Seurat::RunUMAP()` using candidate `umap_neighbors` and `umap_min_dist`

### Fitness Calculation

Fitness is evaluated directly on two-dimensional UMAP coordinates (umap_1, umap_2) combined with cluster identity labels (ident).
The ratio of between-cluster dispersion to within-cluster dispersion is evaluated via the Calinski-Harabasz index using clusterCrit::intCriteria().

### Evolutionary Search Logic

#### Representation
Individuals are encoded as continuous floating-point vectors. In the adaptive mutation scheme, genome length is `2 * number_of_parameters`:
the first half encodes the clustering parameters, while the second half encodes parameter-specific mutation rates.

#### Adaptive Mutation
Mutation rates update multiplicatively using a log-normal random factor. 
Parameter values are mutated with Gaussian noise, scaled by their corresponding mutation rate, and clamped to valid parameter bounds.

#### Selection
The evolutionary algorithm uses tournament selection for parent mating and greedy replacement for survival
under either the **plus** (mu + lambda) or **comma** (mu, lambda) survival strategy.

#### Parallelization
The evolutionary algorithm uses `parallelMap`, at level `ecr.evaluateFitness`, which parallelizes the fitness evaluation for all individuals of the population.
Therefore, each individual can be evaluated in a different parallelization unit.
This is necessary because evaluating one individual means running a full Seurat pipeline, which is relatively expensive.

geneticUMAP uses socket workers for parallelization to enable parallel computation also under Windows operating systems.

The number of processes can be set with the `num_processes` function parameter of `run_optimization()`.

## Important Assumptions & Considerations

### Evaluation Space
Internal clustering metrics are evaluated on 2D UMAP coordinates rather than high-dimensional PCA space or raw expression matrices.
This aligns optimization directly with visual cluster separation in 2D space, though results can be influenced by UMAP embedding distortions.

### Continuous Search Space
Discrete integer parameters (such as PCA dimensions and neighbor counts) are optimized continuously as floating-point values and 
rounded when passed to Seurat routines.

## References
*We retrieved the testdataset included in this package from the Seurat pipeline (https://satijalab.org/seurat/articles/pbmc3k_tutorial) and converted it to h5. 
The original source of the dataset is:
10x Genomics. 3k PBMCs from a Healthy Donor, Single Cell Gene Expression Dataset by Cell Ranger v1.1.0. 2016. url: https://www.10xgenomics.com/welcome?closeUrl=%2Fdatasets&lastTouchOfferName=3k%20PBMCs%20from%20a%20Healthy%20Donor&lastTouchOfferType=Dataset&product=atera&redirectUrl=%2Fdatasets%2F3-k-pbm-cs-from-a-healthy-donor-1-standard-1-1-0