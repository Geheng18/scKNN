#' Example single-cell RNA-seq dataset
#'
#' A dataset containing single-cell RNA-seq data with annotation,
#' count matrix, and phenotype information.
#'
#' @format This dataset contains three objects:
#' \describe{
#'   \item{annotation}{Data frame with 5000 rows (cells) containing cell.type and donor columns}
#'   \item{counts.mat}{Matrix with 5 genes (rows) x 5000 cells (columns)}
#'   \item{pheno.donor}{Named numeric vector of phenotypes for 100 donors}
#' }
#' @source Example data
#' @examples
#' data(scKNN_data)
#' dim(annotation)
#' dim(counts.mat)
#' length(pheno.donor)
"annotation"

#' @rdname annotation
"counts.mat"

#' @rdname annotation
"pheno.donor"
