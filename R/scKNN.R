#' Main Analysis Function for Single-Cell KNN Variance Component Testing
#'
#' @description
#' This is the main function for the scKNN package that performs variance
#' component testing for single-cell RNA-seq data. It automatically tests
#' all cell type and gene region combinations, providing both individual
#' and overall hypothesis tests.
#'
#' @param pheno.donor Numeric vector of phenotype values for each donor
#' @param counts.mat Numeric matrix of gene expression counts (genes x cells)
#' @param annotation Data frame with cell-level annotations containing:
#'   \describe{
#'     \item{cell.type}{Factor indicating the cell type of each cell}
#'     \item{donor}{Factor indicating the donor/individual of each cell}
#'   }
#' @param gene.regions Optional list of gene regions for stratified analysis.
#'   Each element should be a character vector of gene names. If NULL (default),
#'   all genes are treated as one region.
#' @param input.kernel Character vector specifying kernel types (default "product")
#' @param method Character string, either "LMM" or "KNN" (default)
#' @param MINQUE.type Character string, either "MINQUE0" (default) or "MINQUE1"
#' @param lambda.seq Numeric vector of lambda values for cross-validation. If NULL,
#'   cross-validation is skipped and fixed.lambda is used.
#' @param fixed.lambda Numeric value for lambda when skipping cross-validation. Only
#'   used when lambda.seq is NULL. Default is 0 (no regularization).
#' @param n.folds Integer number of CV folds (default = 5)
#' @param cv.criteria CV selection criterion: "MSE" (default), "correlation", or "R2"
#' @param individual.test.type Character string for individual tests:
#'   \describe{
#'     \item{"KNN.ind.test"}{Use KNN.ind.test function (default)}
#'     \item{"KNN.overall.test"}{Use KNN.overall.test function for each component}
#'   }
#' @param design.mat Optional design matrix for fixed effects (covariates)
#' @param verbose Logical, whether to print progress messages (default TRUE)
#' @param significance.level Numeric, significance level for tests (default 0.05)
#'
#' @return A list of class "scKNN.result" containing:
#'   \describe{
#'     \item{individual.pvalues}{Matrix of p-values for individual tests}
#'     \item{overall.pvalue}{P-value for the overall test}
#'     \item{significant.components}{List of significant components}
#'     \item{variance.estimates}{Estimated variance components}
#'     \item{component.info}{Data frame with component descriptions}
#'     \item{optimal.lambda}{Selected lambda value from cross-validation}
#'     \item{cv.results}{Full cross-validation results}
#'     \item{method}{Method used (LMM or KNN)}
#'     \item{n.components.tested}{Number of components tested}
#'   }
#'
#' @examples
#' \dontrun{
#' # Load example data
#' data(scKNN_data)
#'
#' # Basic usage - tests all cell types and all genes
#' result1 <- scKNN(pheno.donor, counts.mat, annotation)
#'
#' # With gene regions
#' regions <- list(pathway1 = c("Gene1", "Gene2", "Gene3"),
#'                 pathway2 = c("Gene4", "Gene5"))
#' result2 <- scKNN(pheno.donor, counts.mat, annotation,
#'                  gene.regions = regions)
#'
#' # Using fixed lambda (no cross-validation)
#' result3 <- scKNN(pheno.donor, counts.mat, annotation,
#'                  lambda.seq = NULL,
#'                  fixed.lambda = 0.1)
#'
#' # Using different individual test type
#' result4 <- scKNN(pheno.donor, counts.mat, annotation,
#'                  individual.test.type = "KNN.overall.test")
#'
#' # View results
#' print(result3)
#' summary(result3)
#' }
#' @importFrom graphics abline barplot legend
#' @export
scKNN <- function(pheno.donor,
                  counts.mat,
                  annotation,
                  gene.regions = NULL,
                  input.kernel = NULL,
                  method = "KNN",
                  MINQUE.type = "MINQUE0",
                  lambda.seq = c(0, 0.01, 0.1, 1),
                  fixed.lambda = 0,
                  n.folds = 5,
                  cv.criteria = "MSE",
                  individual.test.type = c("KNN.ind.test", "KNN.overall.test"),
                  design.mat = NULL,
                  verbose = TRUE,
                  significance.level = 0.05) {

  # Match arguments
  individual.test.type <- match.arg(individual.test.type)

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
  if (!all(c("cell.type", "donor") %in% names(annotation))) {
    stop("annotation must contain 'cell.type' and 'donor' columns")
  }

  # Ensure factors
  if (!is.factor(annotation$cell.type)) {
    annotation$cell.type <- as.factor(annotation$cell.type)
  }
  if (!is.factor(annotation$donor)) {
    annotation$donor <- as.factor(annotation$donor)
  }

  # Get dimensions
  all.cell.types <- levels(annotation$cell.type)
  n.celltypes <- length(all.cell.types)
  n.donors <- length(unique(annotation$donor))
  n.cells <- nrow(annotation)
  n.genes <- nrow(counts.mat)

  # Determine number of regions
  if (is.null(gene.regions)) {
    n.regions <- 1
    region.names <- "all_genes"
  } else {
    n.regions <- length(gene.regions)
    region.names <- names(gene.regions)
    if (is.null(region.names)) {
      region.names <- paste0("region", 1:n.regions)
    }
  }

  if (verbose) {
    cat("scKNN Analysis Pipeline\n")
    cat("======================\n")
    cat("Number of donors:", n.donors, "\n")
    cat("Number of cells:", n.cells, "\n")
    cat("Number of genes:", n.genes, "\n")
    cat("Number of cell types:", n.celltypes, "\n")
    cat("Number of gene regions:", n.regions, "\n")
    cat("Method:", method, "\n")
    cat("Individual test type:", individual.test.type, "\n\n")
  }

  # Step 1: Cross-validation for lambda selection or use fixed lambda
  if (!is.null(lambda.seq)) {
    if (verbose) cat("Step 1: Cross-validation for lambda selection...\n")

    cv.result <- cv.MINQUE(pheno.donor = pheno.donor,
                           counts.mat = counts.mat,
                           annotation = annotation,
                           gene.regions = gene.regions,
                           lambda.seq = lambda.seq,
                           n.folds = n.folds,
                           input.kernel = input.kernel,
                           MINQUE.type = MINQUE.type,
                           method = method,
                           design.mat = design.mat,
                           criteria = cv.criteria,
                           verbose = verbose)

    optimal.lambda <- cv.result$lambda.min

    if (verbose) {
      cat("\nOptimal lambda:", optimal.lambda, "\n\n")
    }
  } else {
    if (verbose) {
      cat("Step 1: Using fixed lambda...\n")
      cat("Fixed lambda:", fixed.lambda, "\n\n")
    }
    optimal.lambda <- fixed.lambda
    cv.result <- NULL
  }

  # Step 2: Fit full model with optimal lambda
  if (verbose) cat("Step 2: Fitting full model with optimal lambda...\n")

  full.model <- MINQUE(pheno.donor = pheno.donor,
                       counts.mat = counts.mat,
                       annotation = annotation,
                       gene.regions = gene.regions,
                       input.kernel = input.kernel,
                       MINQUE.type = MINQUE.type,
                       method = method,
                       lambda = optimal.lambda,
                       design.mat = design.mat)

  # Step 3: Determine component structure
  n.first.order <- n.celltypes * n.regions

  if (method == "KNN") {
    # Components: I (1), J (1), first-order (n.first.order),
    # second-order (n.first.order), interactions (choose(n.first.order, 2))
    n.total.components <- 2 + n.first.order + n.first.order + choose(n.first.order, 2)

    # Indices to test (exclude I and J)
    test.indices <- 3:n.total.components

    # Build component info
    component.info <- data.frame(
      index = 1:n.total.components,
      type = c("Identity", "Intercept",
               rep("First-order", n.first.order),
               rep("Second-order", n.first.order),
               rep("Interaction", choose(n.first.order, 2))),
      name = rep('NULL', n.total.components),
      stringsAsFactors = FALSE
    )

    # Add names for first-order components
    first.order.names <- character(n.first.order)
    idx <- 1
    for (ct in all.cell.types) {
      for (r in 1:n.regions) {
        first.order.names[idx] <- paste0(ct, "_", region.names[r])
        idx <- idx + 1
      }
    }

    # Build all component names
    component.info$name[1:(2 + n.first.order)] <- c("Identity", "Intercept", first.order.names)
    # Add second-order names
    component.info$name[(3 + n.first.order):(2 + 2*n.first.order)] <-
      paste0(first.order.names, "^2")

    # Add interaction names
    if (n.first.order > 1) {
      int.idx <- 3 + 2*n.first.order
      for (i in 1:(n.first.order-1)) {
        for (j in (i+1):n.first.order) {
          component.info$name[int.idx] <- paste0(first.order.names[i],
                                                 " X ",
                                                 first.order.names[j])
          int.idx <- int.idx + 1
        }
      }
    }

  } else { # LMM
    n.total.components <- 1 + n.first.order
    test.indices <- 2:n.total.components

    component.info <- data.frame(
      index = 1:n.total.components,
      type = c("Identity", rep("First-order", n.first.order)),
      name = rep('NULL', n.total.components),
      stringsAsFactors = FALSE
    )

    # Build component names
    first.order.names <- character(n.first.order)
    idx <- 1
    for (ct in all.cell.types) {
      for (r in 1:n.regions) {
        first.order.names[idx] <- paste0(ct, "_", region.names[r])
        idx <- idx + 1
      }
    }

    component.info$name <- c("Identity", first.order.names)
  }

  # Add variance estimates to component info
  component.info$variance.estimate <- full.model$theta

  # Step 4: Individual tests for each component
  if (verbose) cat("\nStep 3: Performing individual tests for each component...\n")

  n.tests <- length(test.indices)
  individual.pvalues <- numeric(n.tests)
  names(individual.pvalues) <- component.info$name[test.indices]

  for (i in 1:n.tests) {
    idx <- test.indices[i]

    if (verbose && i %% 10 == 1) {
      cat("  Testing components", i, "to", min(i+9, n.tests), "of", n.tests, "...\n")
    }

    if (individual.test.type == "KNN.ind.test") {
      test.result <- KNN.ind.test(pheno.donor = pheno.donor,
                                  counts.mat = counts.mat,
                                  annotation = annotation,
                                  test.ind = idx,
                                  gene.regions = gene.regions,
                                  input.kernel = input.kernel,
                                  MINQUE.type = MINQUE.type,
                                  method = method,
                                  lambda = optimal.lambda,
                                  design.mat = design.mat)
    } else {
      # Use KNN.overall.test for single component
      test.result <- KNN.overall.test(pheno.donor = pheno.donor,
                                      counts.mat = counts.mat,
                                      annotation = annotation,
                                      test.ind = idx,
                                      gene.regions = gene.regions,
                                      input.kernel = input.kernel,
                                      MINQUE.type = MINQUE.type,
                                      method = method,
                                      lambda = optimal.lambda,
                                      design.mat = design.mat)
    }

    individual.pvalues[i] <- test.result$pvalue

    # Print significant results as they are found
    if (verbose && test.result$pvalue < significance.level) {
      cat(sprintf("  Found significant: %-30s p = %e",
                  component.info$name[idx],
                  test.result$pvalue))
      if (test.result$pvalue < significance.level/n.tests) {
        cat(" **")
      }
      if (test.result$pvalue < 0.001) {
        cat(" ***")
      }
      cat("\n")
    }
  }

  # Step 5: Overall test
  if (verbose) cat("\nStep 4: Performing overall test for all components...\n")

  overall.test <- KNN.overall.test(pheno.donor = pheno.donor,
                                   counts.mat = counts.mat,
                                   annotation = annotation,
                                   test.ind = test.indices,
                                   gene.regions = gene.regions,
                                   input.kernel = input.kernel,
                                   MINQUE.type = MINQUE.type,
                                   method = method,
                                   lambda = optimal.lambda,
                                   design.mat = design.mat)

  overall.pvalue <- overall.test$pvalue

  # Step 6: Identify significant components
  significant.components <- list(
    uncorrected = names(individual.pvalues)[individual.pvalues < significance.level],
    bonferroni = names(individual.pvalues)[individual.pvalues < significance.level/n.tests],
    fdr = names(individual.pvalues)[p.adjust(individual.pvalues, method = "fdr") < significance.level]
  )

  # Add p-values and significance to component info for ALL components
  component.info$p.value <- NA
  component.info$p.value[test.indices] <- individual.pvalues
  component.info$significant <- ""
  component.info$significant[test.indices] <- ifelse(individual.pvalues < significance.level, "*", "")
  component.info$bonferroni.sig <- ""
  component.info$bonferroni.sig[test.indices] <- ifelse(individual.pvalues < significance.level/n.tests, "*", "")
  component.info$fdr.sig <- ""
  component.info$fdr.sig[test.indices] <- ifelse(p.adjust(individual.pvalues, method = "fdr") < significance.level, "*", "")

  # Create results summary table
  results.table <- data.frame(
    Component = component.info$name[test.indices],
    Type = component.info$type[test.indices],
    Variance = component.info$variance.estimate[test.indices],
    P.value = individual.pvalues,
    Significant = ifelse(individual.pvalues < significance.level, "*", ""),
    Bonferroni = ifelse(individual.pvalues < significance.level/n.tests, "*", ""),
    FDR = ifelse(p.adjust(individual.pvalues, method = "fdr") < significance.level, "*", ""),
    stringsAsFactors = FALSE
  )

  if (verbose) {
    cat("\n======================\n")
    cat("ANALYSIS COMPLETE\n")
    cat("======================\n\n")

    # Print variance component structure with all information
    cat("Variance component structure:\n")
    comp.table <- data.frame(
      index = component.info$index,
      type = component.info$type,
      name = component.info$name,
      variance.estimate = sprintf("%.9f", component.info$variance.estimate),
      p.value = ifelse(is.na(component.info$p.value),
                       "NA",
                       sprintf("%.3e", component.info$p.value)),
      significant = component.info$significant,
      bonferroni = component.info$bonferroni.sig,
      fdr = component.info$fdr.sig,
      stringsAsFactors = FALSE
    )
    print(comp.table, row.names = FALSE)
    cat("\n")

    # Overall test result
    cat("OVERALL TEST\n")
    cat("------------\n")
    cat("P-value:", formatC(overall.pvalue, format = "e", digits = 3))
    if (overall.pvalue < significance.level) {
      cat(" ***SIGNIFICANT***")
    }
    cat("\n\n")

    # Individual tests summary
    cat("INDIVIDUAL TESTS SUMMARY\n")
    cat("------------------------\n")
    cat("Total components tested:", n.tests, "\n")
    cat("Significance threshold:", significance.level, "\n")
    cat("Bonferroni threshold:", formatC(significance.level/n.tests, format = "e", digits = 3), "\n")
    cat("Significant (uncorrected):", length(significant.components$uncorrected),
        sprintf("(%.1f%%)", 100 * length(significant.components$uncorrected) / n.tests), "\n")
    cat("Significant (Bonferroni):", length(significant.components$bonferroni), "\n")
    cat("Significant (FDR):", length(significant.components$fdr), "\n")

    # Component breakdown by type
    cat("\nCOMPONENT BREAKDOWN\n")
    cat("-------------------\n")
    type.table <- table(component.info$type[test.indices])
    for (i in 1:length(type.table)) {
      cat(names(type.table)[i], ":", type.table[i], "components\n")
    }

    # Show top significant results if any
    if (length(significant.components$uncorrected) > 0) {
      cat("\nTOP SIGNIFICANT COMPONENTS\n")
      cat("--------------------------\n")

      sig.idx <- which(individual.pvalues < significance.level)
      sig.pvals <- sort(individual.pvalues[sig.idx])
      n.show <- min(10, length(sig.pvals))

      for (i in 1:n.show) {
        cat(sprintf("%-30s p = %e", names(sig.pvals)[i], sig.pvals[i]))

        # Add significance markers
        if (sig.pvals[i] < significance.level/n.tests) {
          cat(" **")
        }
        if (sig.pvals[i] < 0.001) {
          cat(" ***")
        }
        cat("\n")
      }

      if (length(sig.pvals) > 10) {
        cat("... and", length(sig.pvals) - 10, "more significant components\n")
      }

      cat("\nSignificance codes: *** p<0.001, ** Bonferroni significant\n")
    } else {
      cat("\nNo significant components found at p <", significance.level, "\n")
    }
  }

  # Return results
  result <- list(
    individual.pvalues = individual.pvalues,
    overall.pvalue = overall.pvalue,
    significant.components = significant.components,
    variance.estimates = full.model$theta,
    component.info = component.info,
    results.table = results.table,
    optimal.lambda = optimal.lambda,
    cv.results = cv.result,
    method = method,
    n.components.tested = n.tests,
    test.indices = test.indices,
    n.donors = n.donors,
    n.cells = n.cells,
    n.genes = n.genes,
    n.celltypes = n.celltypes,
    n.regions = n.regions,
    significance.level = significance.level
  )

  class(result) <- "scKNN.result"

  return(result)
}
#' Plot Method for scKNN Results
#'
#' @param x An object of class "scKNN.result"
#' @param type Character, type of plot: "manhattan", "qq", or "variance"
#' @param ... Additional plotting arguments
#'
#' @export
plot.scKNN.result <- function(x, type = c("manhattan", "qq", "variance"), ...) {
  type <- match.arg(type)

  if (type == "manhattan") {
    # Manhattan plot
    plot(-log10(x$individual.pvalues),
         main = "Manhattan Plot of Component P-values",
         xlab = "Component Index",
         ylab = "-log10(p-value)",
         pch = 16,
         col = ifelse(x$individual.pvalues < x$significance.level, "red", "black"),
         ...)

    # Add significance lines
    abline(h = -log10(x$significance.level), col = "blue", lty = 2)
    abline(h = -log10(x$significance.level/x$n.components.tested), col = "red", lty = 2)

    legend("topright",
           c("Not significant", "Significant",
             paste("p =", x$significance.level), "Bonferroni"),
           col = c("black", "red", "blue", "red"),
           pch = c(16, 16, NA, NA),
           lty = c(NA, NA, 2, 2))

  } else if (type == "qq") {
    # Q-Q plot
    expected <- -log10(ppoints(length(x$individual.pvalues)))
    observed <- -log10(sort(x$individual.pvalues))

    plot(expected, observed,
         main = "Q-Q Plot of P-values",
         xlab = "Expected -log10(p)",
         ylab = "Observed -log10(p)",
         pch = 16,
         ...)

    abline(0, 1, col = "red", lty = 2)

  } else if (type == "variance") {
    # Variance components plot
    test.components <- x$component.info[x$test.indices, ]

    barplot(test.components$variance.estimate,
            names.arg = test.components$name,
            las = 2,
            main = "Estimated Variance Components",
            ylab = "Variance",
            col = ifelse(x$individual.pvalues < x$significance.level,
                         "lightblue", "gray"),
            ...)
  }
}
