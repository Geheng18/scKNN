#' Cross-Validation for MINQUE Lambda Selection (Updated)
#'
#' @description
#' Performs k-fold cross-validation to select the optimal regularization parameter
#' (lambda) for MINQUE variance component estimation in single-cell data analysis.
#'
#' @inheritParams MINQUE
#' @param lambda.seq Numeric vector of lambda values to test
#' @param n.folds Integer number of cross-validation folds (default = 5)
#' @param criteria Selection criterion: "MSE", "correlation", or "R2"
#' @param verbose Logical, whether to print progress messages
#'
#' @return A list of class "cv.MINQUE" containing cross-validation results
#' @importFrom stats sd
#' @export
cv.MINQUE <- function(pheno.donor, counts.mat, annotation, gene.regions = NULL,
                      lambda.seq = c(0, 0.01, 0.1, 1),
                      n.folds = 5,
                      input.kernel = NULL,
                      MINQUE.type = 'MINQUE0',
                      method = 'KNN',
                      design.mat = NULL,
                      criteria = c("MSE" ,"correlation", "R2"),
                      verbose = TRUE) {

  # Match criteria argument
  criteria <- match.arg(criteria)

  # Input validation
  if (n.folds < 2) {
    stop("n.folds must be at least 2")
  }
  if (any(lambda.seq < 0)) {
    stop("All lambda values must be non-negative")
  }

  # Get donor information
  donors <- unique(annotation$donor)
  n.donors <- length(donors)

  # Ensure pheno.donor is named with donor IDs
  if (is.null(names(pheno.donor))) {
    if (length(pheno.donor) != n.donors) {
      stop("pheno.donor must have names matching donor IDs or length equal to number of donors")
    }
    names(pheno.donor) <- as.character(donors)
  }

  if (n.donors < n.folds) {
    warning(paste("Number of donors (", n.donors, ") is less than folds (",
                  n.folds, "). Setting folds to ", n.donors))
    n.folds <- n.donors
  }

  # Create fold assignments for donors
  fold.ids <- sample(rep(1:n.folds, length.out = n.donors))
  names(fold.ids) <- donors

  # Initialize storage for results
  n.lambda <- length(lambda.seq)
  fold.errors <- matrix(NA, nrow = n.folds, ncol = n.lambda)
  colnames(fold.errors) <- paste0("lambda_", lambda.seq)

  # Perform cross-validation
  if (verbose) cat("Starting", n.folds, "-fold cross-validation...\n")

  for (fold in 1:n.folds) {
    if (verbose) cat("  Fold", fold, "of", n.folds, "...\n")

    # Split donors into train/test
    test.donors <- donors[fold.ids == fold]
    train.donors <- donors[fold.ids != fold]

    # Get cell indices for train/test sets
    test.cells <- which(annotation$donor %in% test.donors)
    train.cells <- which(annotation$donor %in% train.donors)

    # Split data - ADD drop = FALSE to preserve matrix structure
    counts.train <- counts.mat[, train.cells, drop = FALSE]
    annotation.train <- annotation[train.cells, ]

    # Drop unused levels from all factor columns
    annotation.train <- droplevels(annotation.train)

    # Get phenotypes for train/test
    pheno.train <- pheno.donor[as.character(train.donors)]
    pheno.test <- pheno.donor[as.character(test.donors)]

    # Split design matrix if provided
    if (!is.null(design.mat)) {
      train.idx <- match(train.donors, rownames(design.mat))
      test.idx <- match(test.donors, rownames(design.mat))
      design.train <- design.mat[train.idx, , drop = FALSE]
      design.test <- design.mat[test.idx, , drop = FALSE]
    } else {
      design.train <- NULL
      design.test <- NULL
    }

    # Test each lambda value
    for (i in 1:n.lambda) {
      lambda <- lambda.seq[i]
      # Fit model on training data
      fit <- MINQUE(pheno.train, counts.train, annotation.train, gene.regions = gene.regions,
                    input.kernel = input.kernel,
                    MINQUE.type = MINQUE.type,
                    method = method,
                    lambda = lambda,
                    design.mat = design.train,
                    constrain = FALSE)

      # Predict on test data using combined data
      pred <- predict.MINQUE(y = pheno.train,
                             counts.mat = counts.mat,  # Full data
                             annotation = annotation,   # Full annotation
                             gene.regions = gene.regions,
                             input.kernel = input.kernel,
                             theta = fit$theta,
                             method = method,
                             beta = fit$beta,
                             design.mat.train = design.train,
                             design.mat.test = design.test,
                             train.donors = train.donors,
                             test.donors = test.donors)

      # Calculate error based on criteria
      if (criteria == "MSE") {
        fold.errors[fold, i] <- mean((pheno.test - pred)^2)
      } else if (criteria == "correlation") {
        fold.errors[fold, i] <- -cor(pheno.test, pred)  # Negative for minimization
      } else if (criteria == "R2") {
        ss.res <- sum((pheno.test - pred)^2)
        ss.tot <- sum((pheno.test - mean(pheno.test))^2)
        fold.errors[fold, i] <- -(1 - ss.res/ss.tot)  # Negative for minimization
      }
    }
  }

  # Calculate mean and SE across folds
  cv.errors <- colMeans(fold.errors, na.rm = TRUE)
  cv.se <- apply(fold.errors, 2, sd, na.rm = TRUE) / sqrt(n.folds)
  if (criteria == "MSE") {
    cv.errors.report = cv.errors
  } else if (criteria == "correlation") {
    cv.errors.report = -cv.errors
  } else if (criteria == "R2") {
    cv.errors.report = -cv.errors
  }

  # Find optimal lambda
  index.min <- which.min(cv.errors)
  lambda.min <- lambda.seq[index.min]

  # One-standard-error rule
  threshold <- cv.errors[index.min] + cv.se[index.min]
  index.1se <- max(which(cv.errors <= threshold))
  lambda.1se <- lambda.seq[index.1se]

  # Create result object
  result <- list(
    lambda.seq = lambda.seq,
    cv.errors = cv.errors,
    cv.se = cv.se,
    lambda.min = lambda.min,
    lambda.1se = lambda.1se,
    index.min = index.min,
    index.1se = index.1se,
    fold.errors = fold.errors,
    criteria = criteria,
    n.folds = n.folds
  )

  class(result) <- "cv.MINQUE"

  if (verbose) {
    cat("\nCross-validation complete.\n")
    cat("Using criteria",criteria,"\n")
    cat("Cross-validation error:\n")
    cat(cv.errors.report,'\n')
    cat("Optimal lambda (min):", lambda.min, "\n")
    cat("Optimal lambda (1se):", lambda.1se, "\n")
  }

  return(result)
}
