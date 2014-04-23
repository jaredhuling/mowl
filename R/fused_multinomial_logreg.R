groupMultinomLogReg <- function(x, y, weights = rep(1, nrow(x)), groups = NULL,
                                lambda = NULL, nlambda = 100,
                                intercept = TRUE, verbose = FALSE,
                                irls.maxiter = 30, irls.tol = 1e-10,
                                beta.init = NULL) {
  
  y.f <- as.factor(y)
  if (is.factor(y)) {
    y <- as.numeric(levels(y)[y])
  }
  classes <- levels(y.f)
  K <- length(classes)
  
  nobs <- nrow(x)
  nvars <- ncol(x)
  
  if (is.null(lambda)) {
    invisible(capture.output(suppressWarnings(
      lambda.max <- 2 * lambdamax(x, y, index = groups, model = LinReg(), center = FALSE))))
    
    lambda.min.ratio = ifelse(nobs < nvars, 1e-5, 1e-3)
    
    lambda <- exp( seq(log(as.double(lambda.max)),
                       log(as.double(lambda.min.ratio)),
                       log(as.double(lambda.min.ratio / lambda.max))
                       /nlambda) )
    lambda <- lambda[1:nlambda]
  }
  
  

  len <- if (intercept) {nvars + 1} else {nvars}
  if (!is.null(beta.init)) {
    stopifnot(all(dim(beta.init) == c(K, len)))
  }
  w <- rep(0.5, nobs)
  beta.list <- vector(mode = "list", length = length(lambda))
  betas <- if(is.null(beta.init)) {array(1, dim = c(K, len))} else {beta.init}
  
  for (l in 1:length(lambda)) {
    converged <- rep(FALSE, K)
    for (i in 1:irls.maxiter) {
      prev <- betas
      for (k in 1:K) {
        if (!converged[k]) { 
          y.working <- 1 * (y.f == classes[k]) * sqrt(weights)
          
          init <- if (intercept) {prev[k,-1]} else {prev[k,]}
          
          invisible(capture.output(suppressWarnings(
            beta.tmp <- drop(grplasso(sqrt(weights) * x, y.working, index = groups, center = FALSE, coef.init = init,
                                      model = LinReg(), lambda = lambda[l], weights = w)$coefficients))))
          
          if (intercept) {
            beta[-1] <- beta.tmp
          } else {
            beta <- beta.tmp
          }
          
          
          if (intercept) {
            xwb.tmp <- drop(x %*% beta[-1])
            beta[1] <- mean( y.working - xwb.tmp)
            xwb <- xwb.tmp + beta[1]
          } else {
            xwb <- drop(x %*% beta)
          }
          
          # update weights
          p <- 1 / (1 + exp(-xwb))
          w <- p * (1 - p)
          
          y.working <- xwb + (y - p) / w
          
          betas[k,] <- beta
          
          if (all(abs(beta - prev[k,]) < irls.tol)) {
            converged[k] <- TRUE
          }
          
        }
      }
      if (all(converged)) {
        if (verbose) cat("IRLS Converged at iteration: ", i, "\n")
        break
      }
    }
    beta.list[[l]] <- betas
  }
  beta.list
}



groupFusedMultinomLogReg <- function(x, y, weights, groups = NULL,
                                     lambda = NULL, nlambda = 100,
                                     lambda.lasso = 0, lambda.fused = 0,
                                     intercept = TRUE,
                                     irls.maxiter = 30, irls.tol = 1e-10,
                                     fused.maxiter = 500, fused.tol = 1e-10,
                                     beta.init = NULL) {
  
  y <- as.factor(y)
  
  classes <- levels(y)
  
  gmlr <- groupMultinomLogReg(x, y, weights = weights, groups = groups, 
                              intercept = FALSE, nlambda = 50,
                              irls.maxiter = irls.maxiter, irls.tol = irls.tol)
  
  cat("group lasso stage complete \n")
  
  in.groups <- vector(mode = "list", length = length(gmlr))
  grps <- sort(unique(groups))
  for (i in 1:length(gmlr)) {
    grp.in.tmp <- array(dim = c(nrow(gmlr[[i]]), length(grps)))
    for (g in 1:length(grps)) {
      grp.locs <- which(groups %in% grps[g])
      for (k in 1:nrow(gmlr[[i]])) {
        grp.in <- if(all(abs(gmlr[[i]][k,grp.locs]) < 1e-10)) {
          FALSE
        } else {
          TRUE
        }
        grp.in.tmp[k, g] <- grp.in
      }
    }
    in.groups[[i]] <- grp.in.tmp
  }
  groups.in <- unique(in.groups)
  
  fused.fit <- fusedMultinomLogRegSecond(x, y, weights = weights,
                                         groups = groups, intercept = intercept,
                                         lambda.lasso = lambda.lasso, lambda.fused = lambda.fused,
                                         irls.maxiter = irls.maxiter, irls.tol = irls.tol, 
                                         fused.maxiter = fused.maxiter, fused.tol = fused.tol,
                                         groups.in = groups.in)
  
  structure(list(coefficients = fused.fit, lambda.lasso = lambda.lasso,
                 lambda.fused = lambda.fused, classes = classes), class = "groupSparseFusedFit")
}


fusedMultinomLogRegSecond <- function(x, y, weights = rep(1, nrow(x)), groups = NULL,
                                      lambda.lasso = 0, lambda.fused = 0,
                                      intercept = TRUE,
                                      irls.maxiter = 30, irls.tol = 1e-10, 
                                      fused.maxiter = 500, fused.tol = 1e-10,
                                      beta.init = NULL, groups.in = NULL) {
  
  y.f <- as.factor(y)
  classes <- levels(y.f)
  K <- length(classes)
  G <- length(groups.in)
  
  nobs <- nrow(x)
  nvars <- ncol(x)
  len <- if (intercept) {nvars + 1} else {nvars}
  if (!is.null(beta.init)) {
    stopifnot(all(dim(beta.init) == c(K, len)))
  }
  w <- rep(0.5, nobs)
  betas <- if(is.null(beta.init)) {array(1, dim = c(K, len))} else {beta.init}\
  beta <- betas[1,]
  beta.list <- vector(mode = "list", length = G)
  grps <- sort(unique(groups))
  
  for (g in 1:G) {
    converged <- rep(FALSE, K)
    for (i in 1:irls.maxiter) {
      prev <- betas
      for (k in 1:K) {
        
        if (!converged[k]) { 
          
          which.groups <- grps[which(groups.in[[g]][k, ])]
          in.idx <- which(groups %in% which.groups | is.na(groups))
          groups.current <- groups[in.idx]
          
          y.working <- 1 * (y.f == classes[k]) * sqrt(weights)
          
          init <- if (intercept) {prev[k,-1]} else {prev[k,]}
          
          beta.tmp <- fusedLassoRidge(sqrt(weights) * x[,in.idx], y.working, 
                                      w, groups = groups.current,
                                      lambda.lasso = lambda.lasso, 
                                      lambda.fused = lambda.fused,
                                      maxiter = fused.maxiter,
                                      intercept = FALSE,
                                      tol = fused.tol, beta.init = init[in.idx])
          
          if (intercept) {
            beta[in.idx + 1] <- beta.tmp
          } else {
            beta[in.idx] <- beta.tmp
          }
          
          
          if (intercept) {
            xwb.tmp <- drop(x %*% beta[-1])
            beta[1] <- mean( y.working - xwb.tmp)
            xwb <- xwb.tmp + beta[1]
          } else {
            xwb <- drop(x %*% beta)
          }
          
          # update weights
          p <- 1 / (1 + exp(-xwb))
          w <- p * (1 - p)
          
          y.working <- xwb + (y - p) / w
          
          betas[k,] <- beta
          
          if (all(abs(beta - prev[k,]) < fused.tol)) {
            converged[k] <- TRUE
          }
          
        }
      }
      if (all(converged)) {
        cat("IRLS Converged at iteration: ", i, "\n")
        break
      }
      
    } #end IRLS loop
    beta.list[[g]] <- betas
  }
  beta.list
}

