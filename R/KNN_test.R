#' Overall Hypothesis Test for Variance Components using KNN
#'
#' @description
#' Performs an overall test for multiple variance components in a Kernel Neural Network
#' model using a mixture of chi-square distributions. Tests the null hypothesis that
#' specified variance components are zero.
#'
#' @param pheno.donor Numeric vector of phenotype values for each donor
#' @param counts.mat Numeric matrix of gene expression counts (genes x cells)
#' @param annotation Data frame with cell-level annotations containing cell.type
#'   and donor columns
#' @param gene.regions A list where each element is a character vector of gene names
#' @param test.ind Integer vector of indices specifying which variance components to test
#' @param input.kernel Character vector specifying kernel types for each cell type
#' @param MINQUE.type Character string, either "MINQUE0" or "MINQUE1"
#' @param method Character string, either "LMM" or "KNN" (default)
#' @param lambda Numeric regularization parameter (default = 0)
#' @param weight Optional numeric vector of initial weights for variance components
#' @param design.mat Optional design matrix for fixed effects
#' @param constrain Logical, whether to constrain variance components to be non-negative
#' @param variance.component.list Optional pre-computed variance component matrices
#'
#' @return A list containing:
#'   \describe{
#'     \item{pvalue}{P-value from the mixture of chi-square test}
#'   }
#'
#' @details
#' The test statistic follows a mixture of chi-square distributions under the null
#' hypothesis. The Davies method is used to compute accurate p-values. This test
#' is particularly useful for testing the joint significance of multiple cell types
#' or kernel components.
#'
#' @examples
#' \dontrun{
#' # Test if cell types 2 and 3 contribute to phenotypic variance
#' test_result <- KNN.overall.test(pheno.donor, counts.mat, annotation,
#'                                  test.ind = c(2, 3))
#' print(test_result$pvalue)
#' }
#'
#' @references
#' Davies, R. B. (1980). "Algorithm AS 155: The distribution of a linear
#' combination of chi-squared random variables." Applied Statistics, 323-333.
#'

#' @importFrom CompQuadForm davies
#' @export
KNN.overall.test <- function(pheno.donor, counts.mat, annotation, test.ind, gene.regions = NULL,
                             input.kernel = NULL, MINQUE.type = 'MINQUE0',
                             method = 'KNN', lambda = 0, weight = NULL,
                             design.mat = NULL, constrain = FALSE,
                             variance.component.list = NULL) {

  # Input validation
  if (!is.numeric(test.ind) || any(test.ind < 1)) {
    stop("test.ind must be positive integers")
  }

  # Estimate variance components using MINQUE
  result <- MINQUE(pheno.donor = pheno.donor,
                   counts.mat = counts.mat,
                   annotation = annotation,
                   gene.regions = gene.regions,
                   input.kernel = input.kernel,
                   MINQUE.type = MINQUE.type,
                   method = method,
                   lambda = lambda,
                   weight = weight,
                   design.mat = design.mat,
                   constrain = constrain,
                   variance.component.list = variance.component.list)

  Cmat <- result$Cmat
  theta <- result$theta
  Qmat <- result$Qmat
  variance.component.list <- result$variance.component.list

  # Validate test indices
  if (any(test.ind > length(theta))) {
    stop("test.ind contains indices larger than number of variance components")
  }

  # Compute sigmaH1 (under null hypothesis)
  sigmaH1 <- Reduce(f = "+", x = mapply(FUN = "*", variance.component.list,
                                        Cmat[1,], SIMPLIFY = FALSE))

  # Compute sigmaH2 (under alternative hypothesis)
  sigmaH2.list <- list()
  test.weight <- rep(1, length(test.ind))
  ind <- 1

  for (i in test.ind) {
    sigmaH2 <- Reduce(f = "+", x = mapply(FUN = "*", variance.component.list,
                                          Cmat[i,], SIMPLIFY = FALSE))
    sigmaH2.list[[ind]] <- sigmaH2
    test.weight[ind] <- Cmat[i, i]
    ind <- ind + 1
  }

  # Compute weighted sum
  ratio <- sum(theta[test.ind] * test.weight)
  sigmaH2 <- Reduce(f = "+", x = mapply(FUN = "*", sigmaH2.list, test.weight,
                                        SIMPLIFY = FALSE))
  ratio <- ratio / theta[1]

  # Compute test statistic covariance
  sigma.total <- Qmat %*% (sigmaH2 - ratio * sigmaH1) %*% Qmat

  # Compute eigenvalues and p-value using Davies method
  eigenvalue <- eigen(sigma.total, symmetric = TRUE, only.values = TRUE)$values

  # Remove very small eigenvalues for numerical stability
  # eigenvalue <- eigenvalue[eigenvalue > 1e-15]

  pvalue <- tryCatch({
    davies(q = 0, lambda = eigenvalue)$Qq
  }, error = function(e) {
    # If Liu fails, use Davies
    liu(q = 0, lambda = eigenvalue)
  })

  return(list(pvalue = pvalue))
}


#' Individual Hypothesis Test for a Single Variance Component using KNN
#'
#' @description
#' Tests whether a specific variance component is significantly different from zero
#' using a mixture of chi-square distributions. This is useful for testing the
#' contribution of individual cell types or kernel components.
#'
#' @param pheno.donor Numeric vector of phenotype values for each donor
#' @param counts.mat Numeric matrix of gene expression counts (genes x cells)
#' @param annotation Data frame with cell-level annotations
#' @param test.ind Integer index of the variance component to test
#' @param gene.regions A list where each element is a character vector of gene names
#' @param input.kernel Character vector specifying kernel types
#' @param MINQUE.type Character string, either "MINQUE0" or "MINQUE1"
#' @param method Character string, either "LMM" or "KNN" (default)
#' @param lambda Numeric regularization parameter (default = 0)
#' @param weight Optional numeric vector of initial weights
#' @param design.mat Optional design matrix for fixed effects
#' @param constrain Logical, whether to constrain variance components
#' @param variance.component.list Optional pre-computed variance components
#'
#' @return A list containing:
#'   \describe{
#'     \item{pvalue}{P-value from the mixture of chi-square test}
#'   }
#'
#' @details
#' This function fits both a full model (with all variance components) and a
#' reduced model (without the tested component). The likelihood ratio test
#' statistic follows a mixture of chi-square distributions.
#'
#' @examples
#' \dontrun{
#' # Test if the second variance component (e.g., a cell type) is significant
#' test_result <- KNN.ind.test(pheno.donor, counts.mat, annotation,
#'                             test.ind = 2)
#' print(test_result$pvalue)
#' }
#'
#' @importFrom CompQuadForm davies
#' @export
KNN.ind.test <- function(pheno.donor, counts.mat, annotation, test.ind, gene.regions = NULL,
                         input.kernel = NULL, MINQUE.type = 'MINQUE0',
                         method = 'KNN', lambda = 0, weight = NULL,
                         design.mat = NULL, constrain = FALSE,
                         variance.component.list = NULL) {

  # Input validation
  if (!is.numeric(test.ind) || length(test.ind) != 1 || test.ind < 1) {
    stop("test.ind must be a single positive integer")
  }

  # Estimate variance components under the full model
  result.full <- MINQUE(pheno.donor = pheno.donor,
                        counts.mat = counts.mat,
                        annotation = annotation,
                        gene.regions = gene.regions,
                        input.kernel = input.kernel,
                        MINQUE.type = MINQUE.type,
                        method = method,
                        lambda = lambda,
                        weight = weight,
                        design.mat = design.mat,
                        constrain = constrain,
                        variance.component.list = variance.component.list)

  Cmat <- result.full$Cmat
  theta.full <- result.full$theta
  Qmat <- result.full$Qmat
  variance.component.list.full <- result.full$variance.component.list

  # Validate test index
  if (test.ind > length(theta.full)) {
    stop("test.ind is larger than number of variance components")
  }

  if (test.ind == 1) {
    warning("Testing the residual variance component (index 1) may not be meaningful")
  }

  # Remove the variance component to test to get the reduced model
  variance.component.list.reduce <- variance.component.list.full[-test.ind]

  # Estimate variance components under the reduced model
  result.reduce <- MINQUE(pheno.donor = pheno.donor,
                          counts.mat = counts.mat,
                          annotation = annotation,
                          gene.regions = gene.regions,
                          input.kernel = input.kernel,
                          MINQUE.type = MINQUE.type,
                          method = method,
                          lambda = lambda,
                          weight = weight,
                          design.mat = design.mat,
                          constrain = constrain,
                          variance.component.list = variance.component.list.reduce)

  theta.reduce <- result.reduce$theta

  # Compute sigmaH1 (covariance under reduced model)
  sigmaH1 <- Reduce(f = "+", x = mapply(FUN = "*", variance.component.list.reduce,
                                        theta.reduce, SIMPLIFY = FALSE))

  # Compute sigmaH2 (covariance contribution from tested component)
  sigmaH2 <- Reduce(f = "+", x = mapply(FUN = "*", variance.component.list.full,
                                        Cmat[test.ind,], SIMPLIFY = FALSE))

  # Compute test statistic covariance matrix
  sigmaH1.svd <- svd(sigmaH1)

  # Handle potential numerical issues with small eigenvalues
  tol <- 1e-10
  pos.idx <- which(sigmaH1.svd$d > tol)

  if (length(pos.idx) < length(sigmaH1.svd$d)) {
    warning("sigmaH1 is nearly singular, using pseudoinverse")
    sigmaH1.sqrt <- sigmaH1.svd$u[, pos.idx] %*%
      diag(sqrt(sigmaH1.svd$d[pos.idx])) %*%
      t(sigmaH1.svd$v[, pos.idx])
  } else {
    sigmaH1.sqrt <- sigmaH1.svd$u %*% diag(sqrt(sigmaH1.svd$d)) %*% t(sigmaH1.svd$v)
  }

  sigmaH <- sigmaH1.sqrt %*% Qmat %*% sigmaH2 %*% Qmat %*% sigmaH1.sqrt

  # Compute eigenvalues and p-value using Davies method
  eigenvalue <- eigen(sigmaH, symmetric = TRUE, only.values = TRUE)$values

  # Remove very small eigenvalues for numerical stability
  # eigenvalue <- eigenvalue[eigenvalue > 1e-15]

  pvalue <- tryCatch({
    davies(q = theta.full[test.ind], lambda = eigenvalue)$Qq
  }, error = function(e) {
    # If Liu fails, use Davies
    liu(q = theta.full[test.ind], lambda = eigenvalue)
  })

  return(list(pvalue = pvalue))
}
