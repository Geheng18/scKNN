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
#'   If a single gene is provided (length 1 character vector in a list element),
#'   feature vectors (mean, variance, skewness, etc.) are computed instead of pseudobulk.
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
#' @importFrom stats var median IQR p.adjust ppoints
#' @keywords internal
KNN2LMM <- function(counts.mat, annotation, gene.regions = NULL, input.kernel = NULL, method = 'KNN') {

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
    gene.regions <- list(all.genes)
  } else {
    if (!is.list(gene.regions)) {
      stop("gene.regions must be a list of character vectors")
    }
    all.region.genes <- unlist(gene.regions)
    missing.genes <- setdiff(all.region.genes, all.genes)
    if (length(missing.genes) > 0) {
      stop(paste("The following genes in gene.regions are not in counts.mat:",
                 paste(missing.genes, collapse = ", ")))
    }
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
  pseudobulk.list <- list()
  component.names <- character()

  # Process each cell type and gene region combination
  for (ct in cell.types) {
    for (r in 1:n.regions) {
  
      region.genes <- gene.regions[[r]]
      n.genes.region <- length(region.genes)
  
      # donor x cell indices for this cell type (precompute for speed/consistency)
      ct_cells_idx <- which(annotation$cell.type == ct)
      ann_ct <- annotation[ct_cells_idx, , drop = FALSE]
  
      if (nrow(ann_ct) == 0) {
        # no cells of this type at all -> all zeros
        if (n.genes.region == 1) {
          pseudobulk.mat <- matrix(0, nrow = n.donors, ncol = 6,
                                   dimnames = list(donors, c("mean","var","skew","kurt","median","IQR")))
        } else {
          pseudobulk.mat <- matrix(0, nrow = n.donors, ncol = n.genes.region,
                                   dimnames = list(donors, region.genes))
        }
  
        component.name <- paste0(ct, "_region", r)
        pseudobulk.list[[component.name]] <- pseudobulk.mat
        component.names <- c(component.names, component.name)
        next
      }
  
      # Split indices by donor within this cell type
      donor_split_idx <- split(ct_cells_idx, as.character(annotation$donor[ct_cells_idx]))
  
      # Donor-level library size within (donor, celltype): sum of per-cell lib.size
      lib_sum_by_donor <- vapply(donors, function(d) {
        idx <- donor_split_idx[[d]]
        if (is.null(idx) || length(idx) == 0) return(0)
        sum(annotation$lib.size[idx])
      }, numeric(1))
  
      # -------- CASE 1: Single gene -> features computed on PER-CELL logCPM --------
      if (n.genes.region == 1) {
  
        gene.name <- region.genes[1]
        n.features <- 6
        feature.names <- c("mean", "var", "skew", "kurt", "median", "IQR")
  
        feature.mat <- matrix(0, nrow = n.donors, ncol = n.features,
                              dimnames = list(donors, feature.names))
  
        for (i in 1:n.donors) {
          donor <- donors[i]
          cell.idx <- donor_split_idx[[donor]]
  
          if (!is.null(cell.idx) && length(cell.idx) > 0) {
  
            # raw counts for this gene in this donor x celltype
            gene_counts <- as.numeric(counts.mat[gene.name, cell.idx])
  
            # per-cell logCPM using per-cell lib.size
            lib_cell <- as.numeric(annotation$lib.size[cell.idx])
            # avoid division by 0
            lib_cell <- lib_cell + 1e-8
  
            gene_logcpm <- log2((gene_counts / lib_cell) * 1e6 + 1)
  
            if (length(gene_logcpm) >= 2) {
              feature.mat[i, "mean"]   <- mean(gene_logcpm)
              feature.mat[i, "var"]    <- var(gene_logcpm)
              feature.mat[i, "skew"]   <- skewness(gene_logcpm)
              feature.mat[i, "kurt"]   <- kurtosis(gene_logcpm)
            } else if (length(gene_logcpm) == 1) {
              feature.mat[i, "mean"]   <- gene_logcpm
              feature.mat[i, "var"]    <- 0
              feature.mat[i, "skew"]   <- 0
              feature.mat[i, "kurt"]   <- 0
            }
          }
        }
  
        # standardize features
        feature.mat <- scale(feature.mat, center = TRUE, scale = TRUE)
        feature.mat[is.na(feature.mat)] <- 0
  
        pseudobulk.mat <- feature.mat
  
      } else {
  
        # -------- CASE 2: Multi-gene region -> donor SUM -> donor logCPM --------
        pseudobulk.mat <- matrix(0, nrow = n.donors, ncol = n.genes.region,
                                 dimnames = list(donors, region.genes))
  
        for (i in 1:n.donors) {
          donor <- donors[i]
          cell.idx <- donor_split_idx[[donor]]
  
          if (!is.null(cell.idx) && length(cell.idx) > 0) {
            region.counts <- counts.mat[region.genes, cell.idx, drop = FALSE]
  
            # pseudobulk SUM across cells
            pseudobulk.mat[i, ] <- Matrix::rowSums(region.counts)
          }
        }
  
        # donor-level logCPM using donor-level library size SUM within this celltype
        lib.size <- lib_sum_by_donor + 1e-8
        pseudobulk.mat <- log2(sweep(pseudobulk.mat, 1, lib.size, "/") * 1e6 + 1)
  
        # standardize genes
        pseudobulk.mat <- scale(pseudobulk.mat, center = TRUE, scale = TRUE)
        pseudobulk.mat[is.na(pseudobulk.mat)] <- 0
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

  names(kernel.matrix) <- component.names

  # Construct variance components
  I <- diag(1, n.donors, n.donors)
  J <- matrix(1, n.donors, n.donors)

  variance.component <- c(list(J), kernel.matrix)
  variance.component.list <- list()
  L <- length(kernel.matrix)

  if (method == 'LMM') {
    variance.component.list <- list(Identity = I)
    # Add first-order kernel matrices with names
    for (i in 1:L) {
      variance.component.list[[component.names[i]]] <- kernel.matrix[[i]]
    }
  } else if (method == 'KNN') {
    # Identity and Intercept
    variance.component.list <- list(Identity = I, Intercept = J)

    # First order terms
    for (i in 1:L) {
      variance.component.list[[component.names[i]]] <- kernel.matrix[[i]]
    }

    # Second order terms (squared kernels)
    for (i in 1:L) {
      name <- paste0(component.names[i], "^2")
      variance.component.list[[name]] <- kernel.matrix[[i]]^2
    }

    # Interaction terms
    if (L > 1) {
      for (i in 1:(L-1)) {
        for (j in (i+1):L) {
          name <- paste0(component.names[i], " X ", component.names[j])
          variance.component.list[[name]] <- kernel.matrix[[i]] * kernel.matrix[[j]]
        }
      }
    }
  }

  # Add component information as an attribute for reference
  if (method == 'KNN') {
    # For KNN: I (1), J (1), first-order (L), second-order (L), interactions (choose(L, 2))
    n.total.components <- 2 + L + L + choose(L, 2)

    # Build component info
    component.info <- data.frame(
      index = 1:n.total.components,
      type = c("Identity", "Intercept",
               rep("First-order", L),
               rep("Second-order", L),
               rep("Interaction", choose(L, 2))),
      name = character(n.total.components),
      stringsAsFactors = FALSE
    )

    # Set names for Identity, Intercept, and first-order components
    component.info$name[1:(2 + L)] <- c("Identity", "Intercept", component.names)

    # Add second-order names (squared terms)
    component.info$name[(3 + L):(2 + 2*L)] <- paste0(component.names, "^2")

    # Add interaction names
    if (L > 1) {
      int.idx <- 3 + 2*L
      for (i in 1:(L-1)) {
        for (j in (i+1):L) {
          component.info$name[int.idx] <- paste0(component.names[i], " X ", component.names[j])
          int.idx <- int.idx + 1
        }
      }
    }

    attr(variance.component.list, "component.info") <- component.info

  } else { # LMM
    # For LMM: I (1), first-order (L)
    n.total.components <- 1 + L

    component.info <- data.frame(
      index = 1:n.total.components,
      type = c("Identity", rep("First-order", L)),
      name = c("Identity", component.names),
      stringsAsFactors = FALSE
    )

    attr(variance.component.list, "component.info") <- component.info
  }

  # Keep the other attributes
  attr(variance.component.list, "component.names") <- component.info$name
  attr(variance.component.list, "n.celltypes") <- n.celltypes
  attr(variance.component.list, "n.regions") <- n.regions

  return(variance.component.list)
}

# Helper functions for feature computation
skewness <- function(x) {
  n <- length(x)
  if (n < 3) return(0)
  m <- mean(x)
  s <- sd(x)
  if (s == 0) return(0)
  sum((x - m)^3) / (n * s^3)
}

kurtosis <- function(x) {
  n <- length(x)
  if (n < 4) return(0)
  m <- mean(x)
  s <- sd(x)
  if (s == 0) return(0)
  sum((x - m)^4) / (n * s^4) - 3
}
