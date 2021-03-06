
mowl.fit <- function(x, y, A, groups = NULL, group.sparsity = 0, nfolds, 
                     seed = 123, oracle = NULL, verbose = FALSE, fused = FALSE, lambda = NULL, lambda.min = 1e-2, no.lasso = FALSE,
                     alpha = if(is.null(groups)) {1} else {0.5}, threshold = 0, ...) {
  
  thiscall <- match.call()
  
  nvars <- ncol(x)
  nobs <- nrow(x)
  K <- length(unique(A))
  weights <- y * K
  

  if (is.null(groups)) {
    model <- glmnet(x, A, family = "multinomial", weights = weights, alpha = alpha, lambda = lambda, ...)
  } else {
    config <- msgl.algorithm.config(verbose = FALSE)
    GG <- unique(groups)
    n.groups <- length(GG)
    GG.nonzero <- which(!is.na(GG))
    groups.nonzero <- which(!is.na(groups))
    grouping <- groups; grouping[is.na(grouping)] <- 0
    
    gw <- numeric(length = n.groups); gw[GG.nonzero] <- sqrt(K * (1 - group.sparsity) * (length(GG.nonzero)))
    pw <- array(1, dim = c(K, nvars)); pw[,groups.nonzero] <- 0
    if (no.lasso) pw <- pw * 0
    msgl.lambda <- msgl.lambda.seq(x, A, sampleWeights = weights, groupWeights = gw, grouping = grouping,
                                   parameterWeights = pw, alpha = alpha, d = 100, lambda.min = lambda.min)
    model <- msgl(x, classes = A, sampleWeights = weights, groupWeights = gw, grouping = grouping,
                  parameterWeights = pw, alpha = alpha, lambda = msgl.lambda, algorithm.config = config, ...)
  }
  
  if (threshold > 0) {
    model <- thresholdModel(model, threshold)
  }
  
  if (!is.null(oracle)) {
    
    if (inherits(model, "glmnet")) {
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
  d.vals.cv <- array(0, dim = c(nfolds, length(model$lambda), K))
  
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
                          alpha = alpha, lambda = model$lambda, ...)
    } else {
      fit.fold <- msgl(x.train, classes = A.train, sampleWeights = w.train, groupWeights = gw, 
                       grouping = grouping, parameterWeights = pw, alpha = alpha, 
                       lambda = msgl.lambda, algorithm.config = config, ...)
    }
    
    if (threshold > 0) {
      fit.fold <- thresholdModel(fit.fold, threshold)
    }
    
    if (!is.null(groups)) {
      preds <- predict(fit.fold, x = x.test)$classes
      dimnames(preds) <- NULL
    }
    
    for (i in 1:length(fit.fold$lambda)) {
      if (is.null(groups)) {
        preds <- predict(fit.fold, newx = x.test, type = "class", s = fit.fold$lambda[i])
        values[f, i] <- value.func(A.test, preds, y.test)
        misclass[f, i] <- weighted.mean(preds != A.test, w.test) #mean(preds != oracle.test)
        d.vals.cv[f, i, ] <- computeDfromPreds(preds, y.test, A.test)
      } else {
        preds.i <- drop(preds[,i])
        values[f, i] <- value.func(A.test, preds.i, y.test)
        misclass[f, i] <- weighted.mean(preds.i != A.test, w.test) #mean(preds != oracle.test)
        d.vals.cv[f, i, ] <- computeDfromPreds(preds.i, y.test, A.test)
      }
    }
    if (verbose) cat("Fold = ", f, "\n")
  }
  
  if (!is.null(oracle)) {
    optimal.ind <- which.max(pct.correct)
  }
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
              optimal.d.lambda = model$lambda[optimal.ind.d],
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
              d.aic = d.aic,
              values = values,
              misclass = misclass,
              d.vals = d.vals,
              d.vals.cv = d.vals.cv,
              class.lambda.idx = class.ind,
              value.lambda.idx = value.ind,
              aic.lambda.idx = aic.ind)
  class(ret) <- "owlfit"
  ret
}


computeD <- function(obj, newx, outcome, actual.treatments, group.idx = NULL) {
  
  nlams <- if(inherits(obj, "groupSparseFusedFit2Stage")) {
    length(obj$coefficients[[1]])
  } else if (inherits(obj, "groupSparseFusedFit")) {
    nrow(obj$lambda)
  } else {
    length(obj$lambda)
  }
  
  if (is.factor(actual.treatments)) {
    t.vals <- sort(levels(actual.treatments))
    actual.treatments <- levels(actual.treatments)[actual.treatments]
  } else {t.vals <- sort(unique(actual.treatments))}
  if (is.factor(outcome)) {
    outcome <- levels(outcome)[outcome]
  }
  K <- if (inherits(obj, "msgl")){ 
    obj$beta[[1]]@Dim[1]
  } else if (inherits(obj, "glmnet")){ 
    length(obj$beta)
  } else if (inherits(obj, "groupSparseFusedFit2Stage") | 
               inherits(obj, "groupSparseFusedFit")) {
    length(obj$classes)
  }
  ret <- array(0, dim = c(K, nlams))
  rownames(ret) <- paste("d", 1:K, sep="")
  
  if (inherits(obj, "msgl")) {
    predicted.treatments <- predict(obj, x = newx)$classes
    dimnames(predicted.treatments) <- NULL
  } else if (inherits(obj, "glmnet")){
    predicted.treatments <- predict(obj, newx = newx, s = obj$lambda, type = "class")
  } else if (inherits(obj, "groupSparseFusedFit2Stage")) {
    predicted.treatments <- predict(obj, newx = newx, group.idx = group.idx, type = "class")
  } else if (inherits(obj, "groupSparseFusedFit")) {
    predicted.treatments <- predict(obj, newx = newx, type = "class")
  }
  
  for (i in 1:K) {
    agree.ind <- apply(predicted.treatments, 2, function(x) which(x == actual.treatments & x == t.vals[i]))
    for (l in 1:nlams) {
      ret[i, l] <- if (length(agree.ind) == 0 || length(agree.ind[[l]]) == 0) {NA} else 
                      {mean(outcome[agree.ind[[l]]]) - mean(outcome[-agree.ind[[l]]])}
    }
  }
  ret
}

computeDfromPreds <- function(preds, outcome, actual.treatments) {
  if (is.factor(actual.treatments)) {
    t.vals <- sort(levels(actual.treatments))
    actual.treatments <- levels(actual.treatments)[actual.treatments]
  } else {t.vals <- sort(unique(actual.treatments))}
  if (is.factor(outcome)) {
    outcome <- levels(outcome)[outcome]
  }
  tr.factor <- as.factor(actual.treatments)
  K <- length(levels(tr.factor))

  ret <- numeric(length = K)
  names(ret) <- paste("d", 1:K, sep="")
  names(preds) <- NULL
  dimnames(preds) <- NULL
  
  for (i in 1:K) {
    agree.ind <- which(preds == actual.treatments & preds == t.vals[i])
    ret[i] <- mean(outcome[agree.ind]) - mean(outcome[-agree.ind])
  }
  ret
}

computeDpctCorrect <- function(obj, newx, outcome, actual.treatments, oracle) {
  if (is.factor(actual.treatments)) {
    t.vals <- sort(levels(actual.treatments))
    actual.treatments <- levels(actual.treatments)[actual.treatments]
  } else {t.vals <- sort(unique(actual.treatments))}
  if (is.factor(outcome)) {
    outcome <- levels(outcome)[outcome]
  }
  #K <- if (inherits(obj, "msgl")){obj$beta[[1]]@Dim[1]} else {length(obj$beta)}

  if (inherits(obj, "msgl")) {
    predicted.treatments <- predict(obj, x = newx)$classes
    dimnames(predicted.treatments) <- NULL
    nlams <- length(obj$lambda)
  } else if (inherits(obj, "groupSparseFusedFit")) {
    predicted.treatments <- predict(obj, newx = newx, type = "class")
    nlams <- nrow(obj$lambda)
  } else {
    predicted.treatments <- predict(obj, newx = newx, s = obj$lambda, type = "class")
    nlams <- length(obj$lambda)
  }
  
  a.factor <- as.factor(actual.treatments)
  K <- length(levels(a.factor))
  ret <- array(0, dim = c(K, nlams))
  rownames(ret) <- paste("d", 1:K, sep="")
  
  pct.correct <- apply(predicted.treatments, 2, function(x) mean(x == oracle))
  
  for (i in 1:K) {
    agree.ind <- apply(predicted.treatments, 2, function(x) which(x == actual.treatments & x == t.vals[i]))
    for (l in 1:nlams) {
      ret[i, l] <- mean(outcome[agree.ind[[l]]]) - mean(outcome[-agree.ind[[l]]])
    }
  }
  list(d.vals = ret, pct.correct = pct.correct)
}


computeD.owlfit <- function(obj) {
  K <- if (inherits(obj$model, "msgl")){obj$model$beta[[1]]@Dim[1]} else {length(obj$model$beta)}
  ret <- array(0, dim = c(K, length(obj$model$lambda)))
  rownames(ret) <- paste("d", 1:K)
  for (i in 1:K) {
    if (inherits(obj$model, "msgl")) {
      predict(obj$model, x = x)$classes
    } else {
      predicted.treatments <- predict(obj$model, newx = newx, s = obj$model$lambda, type = "class")
    }
    agree.ind <- apply(predicted.treatments, 2, function(x) which(x == actual.treatments & x == as.character(i)))
    for (l in 1:length(obj$model$lambda)) {
      ret[i, l] <- mean(outcome[agree.ind[[l]]]) - mean(outcome[-agree.ind[[l]]])
    }
  }
  ret
}


