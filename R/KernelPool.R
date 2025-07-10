#' Kernel Matrix Construction Functions
#'
#' @description
#' This file contains a comprehensive collection of kernel functions for constructing
#' similarity matrices from genomic or expression data. Kernels are used to capture
#' different types of relationships between samples and serve as building blocks for
#' variance components in mixed models and kernel methods.
#'
#' @details
#' The kernel functions fall into several categories:
#'
#' **Distance-based kernels:**
#' - `ibs`: Identity-by-state kernel based on allelic similarity
#' - `gaussian`: Gaussian (RBF) kernel using Euclidean distance
#' - `RBF`: Radial basis function kernel
#'
#' **Product-based kernels:**
#' - `product`: Standard inner product kernel
#' - `productidentity`: Product kernel with identity regularization
#' - `absproduct`: Absolute value product kernel
#' - `polynomial`: Polynomial kernel of specified degree
#'
#' **Specialized genomic kernels:**
#' - `car`: Conditional autoregressive kernel accounting for genetic correlation
#' - `delta`: Binary one-hot kernel for categorical data
#'
#' **Transformation kernels:**
#' - `sigmoid`: Sigmoid transformation of product kernel
#' - `relu`: ReLU (rectified linear unit) transformation
#' - `exponential`: Exponential transformation
#' - `SP`: Spherical kernel with angular similarity
#'
#' **Structural kernels:**
#' - `identity`: Identity matrix
#' - `J`: All-ones matrix (intercept kernel)
#' - `sqidentity`: Squared diagonal elements
#' - `fstidentity`: First-order diagonal kernel
#'
#' @name KernelPool
#' @importFrom stats cor dist
NULL

#' Find and Compute Kernel Matrix
#'
#' @description
#' Main dispatcher function that selects and computes the appropriate kernel matrix
#' based on the specified kernel name and parameters.
#'
#' @param kernelName Character string specifying the kernel type. Options include:
#'   "ibs", "car", "identity", "product", "productidentity", "absproduct",
#'   "polynomial", "gaussian", "sigmoid", "relu", "exponential", "J",
#'   "sqidentity", "fstidentity", "RBF", "SP", "delta"
#' @param geno Numeric matrix where rows are samples and columns are features
#'   (e.g., genes, SNPs). For genomic applications, typically coded as 0, 1, 2
#'   for SNP genotypes or expression values for genes.
#' @param maf Numeric vector of minor allele frequencies (used only for "car" kernel)
#' @param weights Numeric vector or scalar for feature weighting (default = 1)
#' @param constant Numeric constant for polynomial kernel (default = 1)
#' @param deg Integer degree for polynomial kernel (default = 2)
#' @param sigma Numeric bandwidth parameter for Gaussian/RBF kernels (default = 1)
#' @param theta Numeric parameter for spherical kernel (default = 1)
#'
#' @return Numeric matrix of size n x n where n is the number of samples (rows in geno).
#'   The matrix is symmetric and represents pairwise similarities between samples.
#'
#' @examples
#' \dontrun{
#' # Example with expression data
#' geno <- matrix(rnorm(100 * 50), nrow = 100, ncol = 50)  # 100 samples, 50 genes
#'
#' # Product kernel
#' K1 <- findKernel("product", geno)
#'
#' # Gaussian kernel with custom bandwidth
#' K2 <- findKernel("gaussian", geno, sigma = 0.5)
#'
#' # Polynomial kernel
#' K3 <- findKernel("polynomial", geno, deg = 3, constant = 1)
#'
#' # Identity kernel
#' K4 <- findKernel("identity", geno)
#' }
#'
#' @seealso
#' Individual kernel functions: \code{\link{getProduct}}, \code{\link{getGaussian}},
#' \code{\link{getIbs}}, etc.
#'
#' @keywords internal
findKernel = function(kernelName, geno, maf = NULL, weights = 1, constant = 1, deg = 2,sigma = 1,theta = 1)
{
  switch(kernelName,
         ibs = {getIbs(geno,weights = 1)},
         car = {getCar(geno = geno, maf = maf, weights = weights)},
         identity = {getIdentity(geno = geno)},
         product = {getProduct(geno = geno, weights = weights)},
         productidentity = {getProductIdentity(geno = geno, weights = weights)},
         absproduct = {getAbsproduct(geno = geno, weights = weights)},
         polynomial = {getPolynomial(geno = geno, weights = weights, constant = constant, deg = deg)},
         gaussian = {getGaussian(geno = geno, weights = weights, sigma = sigma)},
         sigmoid={getSigmoid(geno = geno)},
         relu={getRelu(geno=geno)},
         exponential={getExponential(geno = geno)},
         J = {getJ(geno = geno)},
         sqidentity = {getSqidentity(geno = geno)},
         fstidentity = {getFstidentity(geno = geno)},
         RBF = {getRBF(geno = geno, sigma=sigma)},
         SP = {getSP(geno = geno, theta=theta)},
         delta = { getDelta(geno) },
         stop("Invalied Kernel Name!")
  )
}

#' Delta Kernel for Binary One-Hot Data
#'
#' @description
#' Computes a kernel matrix for binary one-hot encoded data where the kernel
#' value is 1 when two samples share the same "hot" bit (active category).
#'
#' @param geno Binary matrix (0s and 1s) representing one-hot encoded data
#'
#' @return Symmetric binary kernel matrix
#'
#' @keywords internal
getDelta = function(geno) {
  if (!all(geno %in% c(0,1)))
    stop("getDelta: expects a binary one-hot matrix")
  # tcrossprod(geno) gives 1 whenever two rows share the same hot bit
  K = tcrossprod(geno)
  return(K)
}


#' Identity-By-State (IBS) Kernel
#'
#' @description
#' Computes identity-by-state kernel commonly used in genomics. Measures
#' allelic similarity between samples based on shared alleles.
#'
#' @param geno Numeric matrix with genotype data (typically 0, 1, 2 encoding)
#' @param weights Numeric vector of feature weights
#'
#' @return Symmetric IBS kernel matrix
#'
#' @details
#' The IBS kernel counts the number of shared alleles between individuals.
#' For diploid organisms, genotypes are typically coded as 0, 1, 2 representing
#' the number of alternative alleles.
#'
#' @keywords internal
getIbs = function(geno,weights = 1)
{
  n = nrow(geno)
  m = ncol(geno)
  gtemp1 = geno
  gtemp2 = geno
  gtemp1[geno == 2] = 1
  gtemp2[geno == 1] = 0
  gtemp2[geno == 2]= 1
  gtemp = cbind(gtemp1,gtemp2)
  Inner = gtemp %*% diag(weights,nrow = 2*m,ncol = 2*m) %*% t(gtemp)
  X2 = matrix(diag(Inner),nrow = n,ncol = n)
  Y2 = matrix(diag(Inner),nrow = n,ncol = n,byrow = T)
  Bound = sum(matrix(weights,nrow = 1,ncol = m) * 2)
  Dis = (X2 + Y2 - 2 * Inner)
  IBS = Bound - Dis
  return(IBS)
}

#' Conditional Autoregressive (CAR) Kernel
#'
#' @description
#' Computes a conditional autoregressive kernel that accounts for genetic
#' correlation structure. Based on the Cholesky decomposition of a variance
#' matrix that incorporates population structure.
#'
#' @param geno Numeric matrix of genotype or expression data
#' @param maf Numeric vector of minor allele frequencies
#' @param weights Numeric vector of feature weights
#'
#' @return Symmetric CAR kernel matrix
#'
#' @details
#' The CAR kernel is designed to account for population structure and
#' correlation patterns in genomic data. It uses the empirical correlation
#' structure to build a precision matrix.
#'
#' @keywords internal
getCar = function(geno, maf, weights)
{
  S=getIbs(geno,weights)/(2*sum(weights))
  if(weights==1){S=getIbs(geno,weights)/(2*ncol(geno))}
  diag(S)=0
  diag=rowSums(S)
  D=diag(diag)
  gamma = mean(cor(geno))
  Va = D - gamma * S
  VaL = chol(Va)
  K = chol2inv(VaL)
  return(K)
}

#' Identity Kernel
#'
#' @description
#' Returns an identity matrix of appropriate dimensions.
#'
#' @param geno Input matrix (used only for determining dimensions)
#'
#' @return Identity matrix
#'
#' @keywords internal
getIdentity = function(geno)
{
  n = nrow(geno)
  K = diag(rep(1, n))
  return(K)
}

#' Product Kernel (Linear Kernel)
#'
#' @description
#' Computes the standard inner product kernel, also known as linear kernel.
#' Measures linear relationships between samples.
#'
#' @param geno Numeric matrix of features
#' @param weights Numeric vector of feature weights
#'
#' @return Symmetric product kernel matrix
#'
#' @details
#' The product kernel is K(x,y) = <x,y> / (m * sum(weights)) where <x,y>
#' is the inner product and m is the number of features.
#'
#' @keywords internal
getProduct = function(geno, weights = 1)
{
  m = ncol(geno)
  sqrtW = sqrt(diag(weights, m))
  K = (1/(m * sum(weights))) * tcrossprod(geno %*% sqrtW)
  return(K)
}

#' Product Identity Kernel
#'
#' @description
#' Product kernel with added identity regularization to ensure positive definiteness.
#'
#' @param geno Numeric matrix of features
#' @param weights Numeric vector of feature weights
#'
#' @return Regularized product kernel matrix
#'
#' @keywords internal
getProductIdentity = function(geno, weights = 1)
{
  m = ncol(geno)
  n = nrow(geno)
  sqrtW = sqrt(diag(weights, m))
  K = (1/(m * sum(weights))) * tcrossprod(geno %*% sqrtW)
  K = K + 0.2 + diag(0.7,n)
  return(K)
}

#' Absolute Product Kernel
#'
#' @description
#' Product kernel using squared features, emphasizing magnitude relationships.
#'
#' @param geno Numeric matrix of features
#' @param weights Numeric vector of feature weights
#'
#' @return Absolute product kernel matrix
#'
#' @keywords internal
getAbsproduct = function(geno, weights = 1)
{
  m = ncol(geno)
  sqrtW = sqrt(diag(weights, m))
  K = (1/(m * sum(weights))) * tcrossprod((geno^2) %*% sqrtW)
  return(K)
}

#' Polynomial Kernel
#'
#' @description
#' Computes polynomial kernel of specified degree, capturing non-linear
#' polynomial relationships between samples.
#'
#' @param geno Numeric matrix of features
#' @param weights Numeric vector of feature weights
#' @param constant Numeric constant added before exponentiation
#' @param deg Integer degree of the polynomial
#'
#' @return Polynomial kernel matrix
#'
#' @details
#' The polynomial kernel is K(x,y) = (constant + <x,y>)^deg where <x,y>
#' is the normalized inner product.
#'
#' @keywords internal
getPolynomial = function(geno, weights =  1, constant = 1, deg =2)
{
  m = ncol(geno)
  K = (constant + (1/(m * sum(weights))) * tcrossprod(geno)) ^ deg
  return(K)
}

#' Gaussian Kernel (RBF)
#'
#' @description
#' Computes Gaussian (radial basis function) kernel based on Euclidean distances.
#' Captures local similarity patterns with exponential decay.
#'
#' @param geno Numeric matrix of features
#' @param weights Numeric vector of feature weights
#' @param sigma Numeric bandwidth parameter controlling kernel width
#'
#' @return Gaussian kernel matrix
#'
#' @details
#' The Gaussian kernel is K(x,y) = exp(-||x-y||²/(2*m*sigma*sum(weights)))
#' where ||x-y|| is the weighted Euclidean distance.
#'
#' @keywords internal
getGaussian = function(geno, weights = 1, sigma = 1)
{
  m = ncol(geno);
  wtGeno = t(t(geno) * sqrt(weights))
  distMat = as.matrix(dist(wtGeno, method = "euclidean", diag = T))
  K = exp(-(1/(2 * m * (sigma) * sum(weights))) * distMat ^ 2)
  return(K)
}

#' Sigmoid Kernel
#'
#' @description
#' Applies hyperbolic tangent transformation to the product kernel.
#'
#' @param geno Numeric matrix of features
#' @param sigma Numeric shift parameter
#'
#' @return Sigmoid kernel matrix
#'
#' @keywords internal
getSigmoid=function(geno,sigma=1)
{
  K=getProduct(geno)
  K=tanh(K+sigma)
  return(K)
}

#' ReLU Kernel
#'
#' @description
#' Applies rectified linear unit (ReLU) transformation to product kernel.
#' Sets negative values to zero.
#'
#' @param geno Numeric matrix of features
#'
#' @return ReLU kernel matrix
#'
#' @keywords internal
getRelu=function(geno)
{
  K=getProduct(geno)
  K[K<0]=0
  return(K)
}

#' Exponential Kernel
#'
#' @description
#' Applies exponential transformation to the product kernel.
#'
#' @param geno Numeric matrix of features
#'
#' @return Exponential kernel matrix
#'
#' @keywords internal
getExponential=function(geno)
{
  K=getProduct(geno)
  K=exp(K)
  return(K)
}

#' All-Ones Kernel (Intercept Kernel)
#'
#' @description
#' Returns a matrix of all ones, representing a constant intercept term.
#'
#' @param geno Input matrix (used for determining dimensions)
#'
#' @return All-ones matrix
#'
#' @keywords internal
getJ=function(geno)
{
  n = nrow(geno)
  K=matrix(1,n,n)
  return(K)
}

#' Squared Identity Kernel
#'
#' @description
#' Diagonal kernel with squared diagonal elements from the product kernel.
#'
#' @param geno Numeric matrix of features
#'
#' @return Diagonal kernel matrix
#'
#' @keywords internal
getSqidentity=function(geno){
  K = getProduct(geno)
  K = diag(diag(K)^2)
  return(K)
}

#' First-Order Identity Kernel
#'
#' @description
#' Diagonal kernel using diagonal elements from the product kernel.
#'
#' @param geno Numeric matrix of features
#'
#' @return Diagonal kernel matrix
#'
#' @keywords internal
getFstidentity=function(geno){
  K = getProduct(geno)
  K = diag(diag(K))
  return(K)
}

#' Radial Basis Function (RBF) Kernel
#'
#' @description
#' Alternative implementation of RBF kernel using direct Euclidean distances.
#'
#' @param geno Numeric matrix of features
#' @param sigma Numeric bandwidth parameter
#'
#' @return RBF kernel matrix
#'
#' @keywords internal
getRBF=function(geno,sigma=1){
  n = nrow(geno)
  sq_dist = as.matrix(dist(geno,method = 'euclidean'))^2
  K = exp(-sq_dist/(2*sigma^2))
  return(K)
}

#' Spherical Kernel
#'
#' @description
#' Computes spherical kernel based on angular similarity between vectors.
#'
#' @param geno Numeric matrix of features
#' @param theta Numeric parameter controlling sphere radius
#'
#' @return Spherical kernel matrix
#'
#' @details
#' The spherical kernel measures angular similarity and is bounded between 0 and 1.
#'
#' @keywords internal
getSP=function(geno,theta=1){
  num = getProduct(geno)
  row.norm = theta + apply(geno,1,L2.norm)
  deno = sqrt(row.norm%*%t(row.norm))
  K = asin(num/deno)/(2*pi)
  return(K)
}

#' L2 Norm Helper Function
#'
#' @description
#' Computes the L2 norm (mean squared magnitude) of a vector.
#'
#' @param X Numeric vector
#'
#' @return Scalar L2 norm value
#'
#' @keywords internal
L2.norm = function(X){
  m = length(X)
  y = sum(X^2)/m
  return(y)
}
