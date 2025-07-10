#' Predict Method for MINQUE with Single-Cell Data
#'
#' @description
#' Predicts phenotype values for new donors using a fitted MINQUE model. This function
#' computes Best Linear Unbiased Predictions (BLUPs) based on the estimated variance
#' components and observed phenotypes.
#'
#' @param y Numeric vector of training phenotypes (length = number of training donors)
#' @param counts.mat Combined count matrix for both training and test data (genes x cells)
#' @param annotation Combined annotation data frame for both training and test cells
#' @param gene.regions A list where each element is a character vector of gene names
#' @param input.kernel Character vector specifying kernel types for each cell type
#' @param theta Numeric vector of estimated variance components from MINQUE
#' @param method Character string, either "LMM" or "KNN" (must match training)
#' @param beta Numeric vector of fixed effects coefficients (if applicable)
#' @param design.mat.train Design matrix for training donors
#' @param design.mat.test Design matrix for test donors
#' @param variance.component.list Optional pre-computed variance components
#' @param train.donors Character vector of training donor IDs
#' @param test.donors Character vector of test donor IDs
#'
#' @return Numeric vector of predicted phenotype values for test donors
#'
#' @details
#' The prediction uses the mixed model equation:
#' \deqn{\hat{y}_{test} = X_{test}\hat{\beta} + \Sigma_{test,train}\Sigma_{train}^{-1}(y_{train} - X_{train}\hat{\beta})}
#'
#' where \eqn{\Sigma} is the covariance matrix constructed from variance components.
#'
#' @examples
#' \dontrun{
#' # After fitting MINQUE model
#' fit <- MINQUE(pheno.train, counts.train, annotation.train)
#'
#' # Predict for new donors
#' pred <- predict.MINQUE(pheno.train, counts.combined, annotation.combined,
#'                        input.kernel = NULL, theta = fit$theta,
#'                        train.donors = train.ids, test.donors = test.ids)
#' }
#'
#' @importFrom MASS ginv
#' @keywords internal
predict.MINQUE <- function(y, counts.mat, annotation, gene.regions = NULL,
                           input.kernel = NULL, theta, method = 'KNN', beta = NULL,
                           design.mat.train = NULL, design.mat.test = NULL,
                           variance.component.list = NULL,
                           train.donors, test.donors) {

  # Get variance components if not provided
  if (is.null(variance.component.list)) {
    variance.component.list <- KNN2LMM(counts.mat = counts.mat,
                                       annotation = annotation,
                                       gene.regions = gene.regions,
                                       input.kernel = input.kernel,
                                       method = method)
  }

  n.variance.component <- length(variance.component.list)

  # Get dimensions
  N <- nrow(variance.component.list[[1]])  # Total number of donors
  Ntr <- length(train.donors)
  Nts <- length(test.donors)

  # Verify dimensions
  if (length(y) != Ntr) {
    stop("Length of y must equal number of training donors")
  }

  # Initialize fixed effects if not provided
  if (is.null(design.mat.train) || is.null(design.mat.test)) {
    design.mat.train <- matrix(0, Ntr, 1)
    design.mat.test <- matrix(0, Nts, 1)
    beta <- 0
  }

  # Get donor indices in the variance component matrices
  all.donors <- unique(annotation$donor)
  train.idx <- match(train.donors, all.donors)
  test.idx <- match(test.donors, all.donors)

  if (any(is.na(train.idx)) || any(is.na(test.idx))) {
    stop("Some donors not found in variance component matrices")
  }

  # Compute K.U matrix (sum of weighted variance components except residual)
  K.U <- Reduce(`+`, mapply(`*`, variance.component.list[-1], theta[-1], SIMPLIFY = FALSE))

  # Extract relevant submatrices
  if (Nts == 0) {
    # No test data - return fitted values for training data
    sigma <- K.U[train.idx, train.idx]
    Sigma <- K.U[train.idx, train.idx] + diag(theta[1], Ntr)

    y.hat <- as.numeric(design.mat.train %*% beta +
                          sigma %*% ginv(Sigma) %*% (y - design.mat.train %*% beta))
  } else {
    # Predict for test data
    sigma <- K.U[test.idx, train.idx]  # Covariance between test and train
    Sigma <- K.U[train.idx, train.idx] + diag(theta[1], Ntr)  # Training covariance

    y.hat <- as.numeric(design.mat.test %*% beta +
                          sigma %*% ginv(Sigma) %*% (y - design.mat.train %*% beta))
  }

  return(y.hat)
}

