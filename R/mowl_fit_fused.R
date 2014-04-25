

mowl.fit.fused <- function(x, y, A, groups = NULL, group.sparsity = 0, nfolds, 
                           seed = 123, oracle = NULL, verbose = FALSE, ...) {
  
  thiscall <- match.call()
  
  nvars <- ncol(x)
  nobs <- nrow(x)
  K <- length(unique(A))
  weights <- y * K
  

  if (is.null(groups)) {
    stop("You must specify the groups for group, fused lasso")
  }
  model <- groupFusedMultinomLogReg(x, A, groups = groups, weights = weights, ...)
  ngr <- length(model$coefficients)
  nlams <- length(model$coefficients[[1]])
  

  if (!is.null(oracle)) {
    pct.correct <- array(NA, dim = c(ngr, nlams))
    for (g in 1:ngr) {
      preds <- predict(model, newx = x, group.idx = g, type = "class")
      
      pct.correct.tmp <- apply(preds, 2, function(x) mean(x == oracle))
      attr(pct.correct.tmp, "names") <- NULL
      pct.correct[g, ] <- pct.correct.tmp
    }
  }
  
  #if (is.null(groups)) {
  #  aic.fit <- aic(model, y, A, weights)
  #} else {
  #  aic.fit <- aic.msgl(model)
  #}
  #aic.ind <- which.min(aic.fit)
  
  require(dismo)
  set.seed(seed)
  folds <- kfold(y, k = nfolds, by = A)
  misclass <- values <- array(0, dim = c(nfolds, ngr, nlams))
  
  for (f in 1:nfolds) {
    x.train <- x[folds != f,]
    A.train <- A[folds != f]
    w.train <- weights[folds != f]
    x.test <- x[folds == f,]
    A.test <- A[folds == f]
    y.test <- y[folds == f]
    w.test <- weights[folds == f]
    oracle.test <- oracle[folds == f]
    

    fit.fold <- groupFusedMultinomLogReg(x.train, A.train, groups = groups, 
                                         weights = w.train, groups.in = model$groups.in, ...)
    
    
    for (g in 1:ngr) { #loop through group-lasso selections
      
      preds <- predict(fit.fold, newx = x.test, group.idx = g, type = "class")
      for (i in 1:nlams) { #loop through combinations of lasso / fused lasso params

        preds.i <- drop(preds[,i])
        values[f, g, i] <- value.func(A.test, preds.i, y.test)
        misclass[f, g, i] <- weighted.mean(preds.i != A.test, w.test) #mean(preds != oracle.test)
        
      }
    }
    
    if (verbose) cat("Fold = ", f, "\n")
  }
  
  if (!is.null(oracle)) optimal.ind.all <- numeric(ngr)
  class.ind.all <- value.ind.all <- optimal.ind.d.all <- numeric(ngr)
  d.optimal.all <- d.value.all <- d.class.all <- array(NA, dim = c(K, ngr))
  d.vals.all <- vector(mode = "list", length = ngr)
  max.pct.corrects <- numeric(ngr)
  
  for (g in 1:ngr) {
    if (!is.null(oracle)) {
      optimal.ind.tmp <- which.max(pct.correct[g,])
      max.pct.corrects[g] <- pct.correct[g, optimal.ind.tmp]
      optimal.ind.all[g] <- optimal.ind.tmp
    }
    class.ind.tmp <- which.min(colMeans(misclass[,g,]))
    value.ind.tmp <- which.max(colMeans(values[,g,]))
    class.ind.all[g] <- class.ind.tmp
    value.ind.all[g] <- value.ind.tmp
    
    d.vals.tmp <- computeD(model, x, y, A, group.idx = g)
    d.vals.all[[g]] <- d.vals.tmp
    
    optimal.ind.d.tmp <- which.max(colSums(d.vals.tmp))
    optimal.ind.d.all[g] <- if(length(optimal.ind.d.tmp) == 0) {
      NA
    } else {optimal.ind.d.tmp}
    
    d.optimal.tmp <- if(length(optimal.ind.d.tmp) == 0) {
      NA} else {
        d.vals.tmp[, optimal.ind.d.tmp]
      }
    d.value.tmp <- d.vals.tmp[, class.ind.tmp]
    d.class.tmp <- d.vals.tmp[, value.ind.tmp]
    
    d.optimal.all[,g] <- d.optimal.tmp
    d.value.all[,g] <- d.value.tmp
    d.class.all[,g] <- d.class.tmp
    #d.aic <- d.vals[, aic.ind]
  }
  
  d.optimal <- max(d.optimal.all)
  d.value <- max(d.value.all)
  d.class <- max(d.class.all)
  
  if(!is.null(oracle)) g.best.pct <- which.max(max.pct.corrects)
  g.best.optimal <- which.max(d.optimal.all)
  g.best.value <- which.max(d.value.all)
  g.best.class <- which.max(d.class.all)
  
  if(!is.null(oracle)) optimal.ind <- c(g.best.pct, optimal.ind.all[g.best.pct])
  optimal.ind.d <- c(g.best.optimal, optimal.ind.d.all[g.best.optimal])
  class.ind <- c(g.best.class, class.ind.all[g.best.class])
  value.ind <- c(g.best.value, value.ind.all[g.best.value])
  
  ret <- list(model = model,
              call = thiscall,
              lambda = model$lambda,
              optimal.lambda = if(!is.null(oracle)) {model$lambda[optimal.ind[2]]} else {NULL},
              optimal.d.lambda = model$lambda[optimal.ind.d[2],],
              class.lambda = model$lambda[class.ind[2],],
              value.lambda = model$lambda[value.ind[2],],
              #aic.lambda = model$lambda[aic.ind],
              max.pct.correct = if(!is.null(oracle)) {max(pct.correct)} else {NULL},
              class.pct.correct = if(!is.null(oracle)) {pct.correct[class.ind]} else {NULL},
              value.pct.correct = if(!is.null(oracle)) {pct.correct[value.ind]} else {NULL},
              #aic.pct.correct = if(!is.null(oracle)) {pct.correct[aic.ind]} else {NULL},
              d.optimal = d.optimal,
              d.class = d.class,
              d.value = d.value,
              #d.aic = d.aic,
              d.optimal.all,
              d.class.all,
              d.value.all,
              values = values,
              misclass = misclass,
              d.vals = d.vals.all,
              optimal.lambda.idx = optimal.ind.d,
              class.lambda.idx = class.ind,
              value.lambda.idx = value.ind,
              optimal.lambda.idx.all = optimal.ind.d.all,
              class.lambda.idx.all = class.ind.all,
              value.lambda.idx.all = value.ind.all)
              #aic.lambda.idx = aic.ind)
  class(ret) <- c("owlfit", "owlfit.fused")
  ret
}
