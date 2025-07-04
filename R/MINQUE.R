#' Minimum Norm Quadratic Unbiased Estimation (MINQUE) for Single-Cell Data
#'
#' @description
#' Estimates variance components from single-cell RNA-seq data using MINQUE
#' (Minimum Norm Quadratic Unbiased Estimation). This function integrates
#' single-cell expression data with donor-level phenotypes to estimate the
#' contribution of different cell types to phenotypic variation.
#'
#' @param pheno.donor Numeric vector of phenotype values for each donor. Length
#'   must equal the number of unique donors in \code{annotation}.
#' @param counts.mat Numeric matrix of gene expression counts with genes in rows
#'   and cells in columns. Row names should be gene identifiers.
#' @param annotation Data frame with cell-level annotations containing:
#'   \describe{
#'     \item{cell.type}{Factor indicating the cell type of each cell}
#'     \item{donor}{Factor indicating the donor/individual of each cell}
#'   }
#' @param gene.regions A list where each element is a character vector of gene names
#'   defining a gene region. If NULL (default), all genes are treated as one region.
#'   Example: list(c("Gene1", "Gene2"), c("Gene3", "Gene4", "Gene5"))
#' @param input.kernel Character vector specifying kernel types for each cell type.
#'   Options include "product", "gaussian", "polynomial", etc. If NULL (default),
#'   uses "product" kernel for all cell types.
#' @param MINQUE.type Character string specifying the MINQUE estimator type:
#'   \describe{
#'     \item{"MINQUE0"}{Assumes all variance components are equal (default)}
#'     \item{"MINQUE1"}{Uses initial weights for variance components}
#'   }
#' @param method Character string specifying the variance component model:
#'   \describe{
#'     \item{"LMM"}{Linear mixed model with first-order kernels only}
#'     \item{"KNN"}{Kernel neural network with first and second-order terms (default)}
#'   }
#' @param lambda Numeric regularization parameter (default = 0). Larger values
#'   encourage sparsity in variance components.
#' @param weight Numeric vector of initial weights for variance components. Used
#'   only when \code{MINQUE.type = "MINQUE1"}. If NULL, uses equal weights.
#' @param design.mat Optional design matrix for fixed effects (covariates). Each
#'   column represents a covariate, rows correspond to donors.
#' @param constrain Logical indicating whether to constrain variance components
#'   to be non-negative (default = FALSE).
#' @param variance.component.list Optional pre-computed list of variance component
#'   matrices. If NULL, computed internally using \code{KNN2LMM}.
#'
#' @return A list containing:
#'   \describe{
#'     \item{theta}{Numeric vector of estimated variance components}
#'     \item{beta}{Numeric vector of fixed effects coefficients (if design.mat provided)}
#'     \item{Cmat}{Inverse of the C matrix used in estimation}
#'     \item{Qmat}{Precision matrix Q = V^(-1)}
#'     \item{variance.component.list}{List of variance component matrices used}
#'   }
#'
#' @details
#' MINQUE is a method-of-moments approach for estimating variance components
#' without requiring distributional assumptions. The algorithm:
#' \enumerate{
#'   \item Constructs kernel matrices from single-cell data via pseudobulk aggregation
#'   \item Initializes the V matrix based on MINQUE type
#'   \item Computes quadratic forms for estimation
#'   \item Solves for variance components using generalized least squares
#' }
#'
#' The regularization parameter \code{lambda} adds a ridge penalty to all variance
#' components except the residual variance (first component).
#'
#' @examples
#' \dontrun{
#' # Load example data
#' data(counts.mat)
#' data(annotation)
#' data(pheno.donor)
#'
#' # Basic usage with default settings
#' result <- MINQUE(pheno.donor, counts.mat, annotation)
#' print(result$theta)
#'
#' # With covariates
#' design.mat <- cbind(1, rnorm(length(pheno.donor))) # intercept + covariate
#' result_cov <- MINQUE(pheno.donor, counts.mat, annotation,
#'                      design.mat = design.mat)
#'
#' # Using MINQUE1 with custom weights
#' n_components <- length(unique(annotation$cell.type)) + 1
#' custom_weights <- rep(1/n_components, n_components)
#' result_minque1 <- MINQUE(pheno.donor, counts.mat, annotation,
#'                          MINQUE.type = "MINQUE1",
#'                          weight = custom_weights)
#' }
#' \dontrun{
#' # Advanced example with regularization
#' # Test different lambda values
#' lambdas <- c(0, 0.1, 1, 10)
#' results <- lapply(lambdas, function(l) {
#'   MINQUE(pheno.donor, counts.mat, annotation,
#'          lambda = l, constrain = TRUE)
#' })
#'
#' # Compare variance component estimates
#' theta_mat <- sapply(results, function(x) x$theta)
#' matplot(lambdas, t(theta_mat), type = "b",
#'         xlab = "Lambda", ylab = "Variance Component",
#'         main = "Regularization Path")
#' }
#'
#' @references
#' Rao, C.R. (1971). "Estimation of variance and covariance components - MINQUE
#' theory." Journal of Multivariate Analysis, 1(3), 257-275.
#'
#' Zhu, Z., & Zhang, X. (2020). "Penalized variance component analysis of
#' complex traits with single-cell data." Nature Communications, 11(1), 1-12.
#'
#' @seealso
#' \code{\link{KNN2LMM}} for variance component construction,
#' \code{\link{findKernel}} for kernel computation
#'
#' @importFrom Matrix nearPD
#' @export
MINQUE <- function(pheno.donor, counts.mat, annotation, gene.regions = NULL,
                   input.kernel = NULL,MINQUE.type = 'MINQUE0', method = 'KNN',
                   lambda = 0,weight = NULL, design.mat = NULL, constrain = FALSE,
                   variance.component.list = NULL) {

  # Input validation
  if (!is.numeric(pheno.donor)) {
    stop("pheno.donor must be a numeric vector")
  }
  if (!is.matrix(counts.mat)) {
    stop("counts.mat must be a matrix")
  }
  if (!is.data.frame(annotation)) {
    stop("annotation must be a data frame")
  }
  if (!MINQUE.type %in% c("MINQUE0", "MINQUE1")) {
    stop("MINQUE.type must be either 'MINQUE0' or 'MINQUE1'")
  }
  if (lambda < 0) {
    stop("lambda must be non-negative")
  }

  # Check donor consistency
  n_donors <- length(unique(annotation$donor))
  if (length(pheno.donor) != n_donors) {
    stop(paste("Length of pheno.donor (", length(pheno.donor),
               ") must equal number of unique donors (", n_donors, ")"))
  }

  # If variance components are not provided, compute them
  if (is.null(variance.component.list)) {
    variance.component.list = KNN2LMM(counts.mat = counts.mat,
                                      annotation = annotation,
                                      gene.regions = gene.regions,
                                      input.kernel = input.kernel,
                                      method = method)
  }
  n.variance.component = length(variance.component.list)
  N = length(pheno.donor)

  # Initialize V matrix based on MINQUE type
  if (MINQUE.type == 'MINQUE0') {
    Vmat = diag(1, N, N)
  } else if (MINQUE.type == 'MINQUE1') {
    if (is.null(weight)) {
      weight = rep(1 / n.variance.component, n.variance.component)
    }
    if (length(weight) != n.variance.component) {
      stop(paste("Length of weight must equal number of variance components (", n.variance.component, ")"))
    }
    Vmat = Reduce(`+`, mapply(`*`, variance.component.list, weight, SIMPLIFY = FALSE))
  }

  # Ensure Vmat is positive definite
  NPD = nearPD(Vmat)
  Vmat = as.matrix(NPD$mat)
  Qmat = solve(Vmat)

  # Compute Rv and Ry
  Rv = list()
  if (is.null(design.mat)) {
    # No fixed effects
    for (i in seq.int(n.variance.component)) {
      Rv[[i]] = Qmat %*% variance.component.list[[i]]
    }
    Ry = Qmat %*% pheno.donor
  } else {
    # With fixed effects
    if (nrow(design.mat) != N) {
      stop("Number of rows in design.mat must equal length of pheno.donor")
    }
    Vmat.inv = solve(Vmat)
    Pmat = design.mat %*% solve(t(design.mat) %*% Vmat.inv %*% design.mat) %*% t(Vmat.inv)
    for (i in seq.int(n.variance.component)) {
      Rv[[i]] = Qmat %*% (variance.component.list[[i]] - Pmat %*% variance.component.list[[i]])
    }
    Ry = Qmat %*% (pheno.donor - Pmat %*% pheno.donor)
  }

  # Compute U matrix
  Umat = numeric(n.variance.component)
  for (i in seq.int(n.variance.component)) {
    Umat[i] = sum(Ry * (variance.component.list[[i]] %*% Ry))
  }

  # Compute C matrix
  Cmat = matrix(0, n.variance.component, n.variance.component)
  for (i in seq.int(n.variance.component)) {
    for (j in seq.int(n.variance.component)) {
      Cmat[i, j] = sum(Rv[[i]] * Rv[[j]])
    }
  }

  # Adjust C matrix with penalty
  penalty = diag(lambda, n.variance.component, n.variance.component)
  penalty[1, 1] = 0  # No penalty on the first variance component
  Cmat = Cmat + penalty

  # Ensure Cmat is positive definite
  NPD = nearPD(Cmat)
  Cmat = as.matrix(NPD$mat)

  # Estimate variance components
  Cmat = solve(Cmat)
  theta = Cmat %*% Umat

  # Estimate fixed effects if design.mat is provided
  if (!is.null(design.mat)) {
    beta = solve(t(design.mat) %*% Qmat %*% design.mat) %*% t(design.mat) %*% Qmat %*% pheno.donor
    if (!is.null(colnames(design.mat))) {
      names(beta) = colnames(design.mat)
    }
  } else {
    beta = NULL
  }

  theta = as.numeric(theta)
  if (constrain) {
    theta[theta < 0] = 0
  }

  return(list(theta = theta, beta = beta, Cmat = Cmat, Qmat = Qmat,
              variance.component.list = variance.component.list))
}
