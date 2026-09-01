# Backward purging for Random Forest classification.
# At each step, repeated forests stabilize permutation importance before the
# two least-important predictors are removed. Rows are indexed by the number
# of predictors remaining, so a requested model size can be selected by its
# actual variable count rather than by an implicit row offset.

mcRFop_cls <- function(x, y, nTree = 100L, nRep = 10L, seed = NULL) {
  x <- as.data.frame(x)
  y <- factor(y)

  if (ncol(x) < 3L) stop("At least three predictors are required.")
  if (nlevels(y) < 2L) stop("The response must contain at least two classes.")

  detected_cores <- parallel::detectCores(logical = FALSE)
  if (is.na(detected_cores)) detected_cores <- 1L
  n_core <- min(
    max(1L, detected_cores - 1L),
    as.integer(nTree)
  )
  ntree_vec <- rep(nTree %/% n_core, n_core)
  ntree_vec[seq_len(nTree %% n_core)] <-
    ntree_vec[seq_len(nTree %% n_core)] + 1L
  ntree_vec <- ntree_vec[ntree_vec > 0L]

  cl <- parallel::makeCluster(length(ntree_vec), type = "SOCK")
  doSNOW::registerDoSNOW(cl)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  if (!is.null(seed)) parallel::clusterSetRNGStream(cl, iseed = seed)

  selected <- data.frame(
    Accy = rep(NA_real_, ncol(x)),
    variable = rep(NA_character_, ncol(x)),
    stringsAsFactors = FALSE
  )
  x_current <- x

  repeat {
    imp_sum <- NULL
    rf <- NULL

    for (r in seq_len(nRep)) {
      rf <- foreach::foreach(
        ntree = ntree_vec,
        .combine = randomForest::combine,
        .packages = "randomForest"
      ) %dopar% {
        randomForest::randomForest(
          x_current,
          y,
          ntree = ntree,
          importance = TRUE
        )
      }

      imp_r <- randomForest::importance(rf)
      score_col <- if ("MeanDecreaseAccuracy" %in% colnames(imp_r)) {
        "MeanDecreaseAccuracy"
      } else {
        colnames(imp_r)[max(1L, ncol(imp_r) - 1L)]
      }
      score <- imp_r[, score_col]
      imp_sum <- if (is.null(imp_sum)) score else imp_sum + score
    }

    mean_imp <- imp_sum / nRep
    ranked <- names(sort(mean_imp, decreasing = FALSE))
    n_var <- length(ranked)
    selected$Accy[n_var] <- mean(predict(rf, x_current) == y)
    selected$variable[n_var] <- toString(ranked)

    if (n_var <= 3L) break
    x_current <- x_current[, !names(x_current) %in% ranked[1:2], drop = FALSE]
  }

  selected
}
