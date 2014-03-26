
mowl.fit <- function(x, y, A, groups = NULL, nfolds, seed = 123, oracle = NULL, verbose = FALSE) {
  
  thiscall <- match.call()
  weights <- y * 3
  nvars <- ncol(x)
  nobs <- nrow(x)
  K <- length(unique(A))
  
  if (is.null(groups)) {
    model <- glmnet(x, A, family = "multinomial", weights = weights, alpha = 1)
  } else {
    gw <- numeric(length = nvars); gw[groups] <- sqrt(K) * (p - length(groups))
    pw <- array(1, dim = c(K, nvars)); pw[,groups] <- 0
    msgl.lambda <- msgl.lambda.seq(x, A, sampleWeights = weights, groupWeights = gw, 
                                   parameterWeights = pw, alpha = .5, d = 100, lambda.min = 1e-3)
    model <- msgl(x, classes = A, sampleWeights = weights, groupWeights = gw, 
                  parameterWeights = pw, alpha = 0.5, lambda = msgl.lambda)
  }
  
  if (!is.null(oracle)) {
    if (is.null(groups)) {
      preds <- predict(model, newx = x, type = "class", s = model$lambda)
    } else {
      preds <- predict(model, x = x)$classes
    }
    pct.correct <- apply(preds, 2, function(x) mean(x == oracle))
    attr(pct.correct, "names") <- NULL
  }
  
  if (is.null(groups)) {
    aic.fit <- aic(model, y, A, weights)
  } else {
    aic.fit <- aic.msgl(model)
  }
  aic.ind <- which.min(aic.fit)
  
  require(dismo)
  set.seed(seed)
  folds <- kfold(y, k = nfolds, by = A)
  misclass <- values <- array(0, dim = c(nfolds, length(model$lambda)))
  
  for (f in 1:nfolds) {
    x.train <- x[folds != f,]
    A.train <- A[folds != f]
    w.train <- weights[folds != f]
    x.test <- x[folds == f,]
    A.test <- A[folds == f]
    y.test <- y[folds == f]
    w.test <- weights[folds == f]
    oracle.test <- oracle[folds == f]
    
    if (is.null(groups)) {
      fit.fold <- glmnet(x.train, A.train, family = "multinomial", weights = w.train, 
                          alpha = 1, lambda = model$lambda)
    } else {
      fit.fold <- msgl(x.train, classes = A.train, sampleWeights = weights, groupWeights = gw, 
                       parameterWeights = pw, alpha = 0.5, lambda = msgl.lambda)
    }
    
    for (i in 1:length(fit.fold$lambda)) {
      if (is.null(groups)) {
        preds <- predict(fit.fold, newx = x.test, type = "class", s = gfit.fold$lambda[i])
      } else {
        preds <- predict(fit.fold, x = x.test)$classes
      }
      values[f, i] <- value.func(A.test, preds, y.test)
      misclass[f, i] <- weighted.mean(preds != A.test, w.test) #mean(preds != oracle.test)
    }
    if (verbose) cat("Fold = ", f, "\n")
  }
  optimal.ind <- which.max(pct.correct)
  class.ind <- which.min(colMeans(misclass))
  value.ind <- which.max(colMeans(values))
  
  
  d.vals <- computeD(model, x, y, A)
  optimal.ind.d <- which.max(colSums(d.vals))
  
  d.optimal <- d.vals[, optimal.ind.d]
  d.value <- d.vals[, value.ind]
  d.class <- d.vals[, class.ind]
  d.aic <- d.vals[, aic.ind]
  
  ret <- list(model = model,
              call = thiscall,
              optimal.lambda = if(!is.null(oracle)) {model$lambda[optimal.ind]} else {NULL},
              optimal.d.lambda = if(!is.null(oracle)) {model$lambda[optimal.ind.d]} else {NULL},
              class.lambda = model$lambda[class.ind],
              value.lambda = model$lambda[value.ind],
              aic.lambda = model$lambda[aic.ind],
              max.pct.correct = if(!is.null(oracle)) {max(pct.correct)} else {NULL},
              class.pct.correct = if(!is.null(oracle)) {pct.correct[class.ind]} else {NULL},
              value.pct.correct = if(!is.null(oracle)) {pct.correct[value.ind]} else {NULL},
              aic.pct.correct = if(!is.null(oracle)) {pct.correct[aic.ind]} else {NULL},
              d.optimal = d.optimal,
              d.class = d.class,
              d.value = d.value,
              d.aic = d.aic)
  class(ret) <- "owlfit"
  ret
}


computeD <- function(obj, newx, outcome, actual.treatments) {
  ret <- array(0, dim = c(length(obj$beta), length(obj$lambda)))
  rownames(ret) <- paste("d", 1:length(obj$beta))
  for (i in 1:length(obj$beta)) {
    predicted.treatments <- predict(obj, newx = newx, s = obj$lambda, type = "class")
    agree.ind <- apply(predicted.treatments, 2, function(x) which(x == actual.treatments & x == as.character(i)))
    for (l in 1:length(obj$lambda)) {
      ret[i, l] <- mean(outcome[agree.ind[[l]]]) - mean(outcome[-agree.ind[[l]]])
    }
  }
  ret
}

computeD.owlfit <- function(obj) {
  ret <- array(0, dim = c(length(obj$model$beta), length(obj$model$lambda)))
  rownames(ret) <- paste("d", 1:length(obj$model$beta))
  for (i in 1:length(obj$model$beta)) {
    predicted.treatments <- predict(obj$model, newx = newx, s = obj$lambda, type = "class")
    agree.ind <- apply(predicted.treatments, 2, function(x) which(x == actual.treatments & x == as.character(i)))
    for (l in 1:length(obj$model$lambda)) {
      ret[i, l] <- mean(outcome[agree.ind[[l]]]) - mean(outcome[-agree.ind[[l]]])
    }
  }
  ret
}


