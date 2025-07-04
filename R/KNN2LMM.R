#' Convert Single-Cell Data to Kernel Variance Components with Gene Regions
#'
#' @description
#' Transforms single-cell RNA-seq data into kernel-based variance components for
#' linear mixed models. This function performs pseudobulk aggregation by cell type
#' and gene region, then constructs kernel matrices that can be used as variance
#' components in downstream association analyses.
#'
#' @param counts.mat A numeric matrix of gene expression counts with genes in rows
#'   and cells in columns. Row names should be gene identifiers.
#' @param annotation A data frame with cell-level annotations containing at least
#'   two columns:
#'   \describe{
#'     \item{cell.type}{Factor indicating the cell type of each cell}
#'     \item{donor}{Factor indicating the donor/individual of each cell}
#'   }
#'   The number of rows must equal the number of columns in \code{counts.mat}.
#' @param gene.regions A list where each element is a character vector of gene names
#'   defining a gene region. If NULL (default), all genes are treated as one region.
#'   Example: list(c("Gene1", "Gene2"), c("Gene3", "Gene4", "Gene5"))
#' @param input.kernel Character vector specifying kernel types. If gene.regions is
#'   provided, length should equal length(gene.regions) * n_celltypes. Otherwise,
#'   length should equal n_celltypes. If NULL, uses "product" kernel for all.
#' @param method Character string specifying the model type:
#'   \describe{
#'     \item{"LMM"}{Linear mixed model with only first-order kernels}
#'     \item{"KNN"}{Kernel neural network with first-order, second-order, and
#'       interaction terms (default)}
#'   }
#'
#' @return A list of variance component matrices. The order of components (after
#'   the identity matrix) follows: celltype1_region1, celltype1_region2, ...,
#'   celltype2_region1, celltype2_region2, etc.
#'
#' @export
KNN2LMM <- function(counts.mat, annotation, gene.regions = NULL,input.kernel = NULL, method = 'KNN') {

  # Input validation
  if (!is.matrix(counts.mat)) {
    stop("counts.mat must be a matrix")
  }
  if (!is.data.frame(annotation)) {
    stop("annotation must be a data frame")
  }
  if (!all(c("cell.type", "donor") %in% names(annotation))) {
    stop("annotation must contain 'cell.type' and 'donor' columns")
  }
  if (ncol(counts.mat) != nrow(annotation)) {
    stop("Number of cells in counts.mat must match rows in annotation")
  }
  if (!method %in% c("LMM", "KNN")) {
    stop("method must be either 'LMM' or 'KNN'")
  }

  # Handle gene regions
  all.genes <- rownames(counts.mat)
  if (is.null(gene.regions)) {
    # If no gene regions specified, use all genes as one region
    gene.regions <- list(all.genes)
  } else {
    # Validate gene regions
    if (!is.list(gene.regions)) {
      stop("gene.regions must be a list of character vectors")
    }

    # Check that all specified genes exist in counts.mat
    all.region.genes <- unlist(gene.regions)
    missing.genes <- setdiff(all.region.genes, all.genes)
    if (length(missing.genes) > 0) {
      stop(paste("The following genes in gene.regions are not in counts.mat:",
                 paste(missing.genes, collapse = ", ")))
    }

    # Check for overlapping genes
    if (length(all.region.genes) != length(unique(all.region.genes))) {
      stop("Gene regions contain overlapping genes. Each gene should appear in only one region.")
    }
  }

  # Get unique cell types and donors
  cell.types <- levels(annotation$cell.type)
  donors <- levels(annotation$donor)
  n.donors <- length(donors)
  n.celltypes <- length(cell.types)
  n.regions <- length(gene.regions)

  # Initialize list to store pseudobulk matrices
  # Will store as: celltype1_region1, celltype1_region2, ..., celltype2_region1, ...
  pseudobulk.list <- list()
  component.names <- character()

  # Process each cell type and gene region combination
  for (ct in cell.types) {
    for (r in 1:n.regions) {
      # Get genes for this region
      region.genes <- gene.regions[[r]]
      n.genes.region <- length(region.genes)

      # Initialize pseudobulk matrix for this cell type and region
      pseudobulk.mat <- matrix(0, nrow = n.donors, ncol = n.genes.region,
                               dimnames = list(donors, region.genes))

      # Process each donor
      for (i in 1:n.donors) {
        donor <- donors[i]
        # Get indices of cells for this donor and cell type
        cell.idx <- which(annotation$donor == donor & annotation$cell.type == ct)

        if (length(cell.idx) > 0) {
          # Extract counts for the specific genes in this region
          region.counts <- counts.mat[region.genes, cell.idx, drop = FALSE]

          # Calculate mean expression across cells
          if (length(cell.idx) == 1) {
            pseudobulk.mat[i,] <- region.counts
          } else {
            pseudobulk.mat[i,] <- rowMeans(region.counts)
          }
        }
      }

      # Log transform and scale
      pseudobulk.mat <- log2(pseudobulk.mat + 1)
      pseudobulk.mat <- scale(pseudobulk.mat, center = TRUE, scale = TRUE)

      # Handle genes with zero variance
      zero.var <- which(is.na(pseudobulk.mat[1,]))
      if (length(zero.var) > 0) {
        pseudobulk.mat[, zero.var] <- 0
      }

      # Add to list with descriptive name
      component.name <- paste0(ct, "_region", r)
      pseudobulk.list[[component.name]] <- pseudobulk.mat
      component.names <- c(component.names, component.name)
    }
  }

  # Generate kernel matrices
  n.components <- length(pseudobulk.list)
  kernel.matrix <- list()

  if (is.null(input.kernel)) {
    input.kernel <- rep('product', n.components)
  } else if (length(input.kernel) == 1) {
    input.kernel <- rep(input.kernel, n.components)
  } else if (length(input.kernel) != n.components) {
    stop(paste("input.kernel must have length 1 or", n.components,
               "(n_celltypes * n_regions)"))
  }

  for (i in 1:n.components) {
    km.temp <- findKernel(input.kernel[i], pseudobulk.list[[i]])
    kernel.matrix <- c(kernel.matrix, list(km.temp))
  }

  # Add names to kernel matrices for clarity
  names(kernel.matrix) <- component.names

  # Construct variance components
  I <- list(diag(1, n.donors, n.donors))
  J <- list(matrix(1, n.donors, n.donors))

  # variance.component = [J, all raw kernels]
  variance.component <- c(J, kernel.matrix)

  # Final list to collect every variance component
  variance.component.list <- list()
  L <- length(kernel.matrix)  # number of "raw" kernels

  if (method == 'LMM') {
    variance.component.list <- c(I, kernel.matrix)
  } else if (method == 'KNN') {
    # First order terms
    variance.component.list <- c(I, variance.component)

    # Second order terms (squared kernels)
    for (i in 2:(L+1)) {
      variance.component.temp <- list(variance.component[[i]]^2)
      variance.component.list <- c(variance.component.list, variance.component.temp)
    }

    # Interaction terms
    if (L > 1) {
      for (i in 2:(L+1)) {
        if (i < L+1) {
          for (j in (i+1):(L+1)) {
            variance.component.temp <- list(variance.component[[i]] * variance.component[[j]])
            variance.component.list <- c(variance.component.list, variance.component.temp)
          }
        }
      }
    }
  }

  # Add component information as an attribute for reference
  attr(variance.component.list, "component.names") <- c("Identity", "Intercept", component.names)
  attr(variance.component.list, "n.celltypes") <- n.celltypes
  attr(variance.component.list, "n.regions") <- n.regions

  return(variance.component.list)
}
