

ridge.fit <- function(x, y, w = rep(1, nrow(x)), lambda, mult.factor = NULL, intercept = TRUE) {
  
  nobs <- nrow(x)
  cwmx <- colWeightedMeans(x, w)
  x.tilde <- x - matrix(rep(cwmx, nobs), ncol=ncol(x), byrow=TRUE)
  y.tilde <- y - weighted.mean(y, w)
  s <- colWeightedMeans(x^2, w)
  
  if (!is.null(mult.factor)) {
    A <- crossprodcpp(x.tilde, drop(w)) + s * sqrt(nobs) * mult.factor
  } else {
    A <- crossprodcpp(x.tilde, drop(w)) + lambda * sqrt(nobs) * diag(s)
  }
  beta <- drop(solve(A, crossprod(x.tilde, drop(w) * y.tilde)))
  if (intercept) {
    beta <- c(mean(y - x %*% beta), beta)
  }
  beta
}

fusedLassoRidge <- function(x, y, w = rep(1, nrow(x)), groups = NULL, lambda.lasso = 0, lambda.fused = 0, 
                            maxiter = 1000, tol = 1e-10, intercept = TRUE, beta.init = NULL,
                            verbose = FALSE) {
  
  beta <- if (is.null(beta.init)) {
    len <- if (intercept) {ncol(x) + 1} else {ncol(x)}
    rep(1, len)
  } else {
    if (intercept) {
      stopifnot(length(beta.init) == ncol(x) + 1)
    } else {
      stopifnot(length(beta.init) == ncol(x))
    }
    beta.init
  }
  nvars <- ncol(x)
  for (i in 1:maxiter) {
    prev <- beta
    beta.pen <- if (intercept) {prev[-1]} else {prev}
    
    A <- lambda.lasso * diag(1 / abs(beta.pen))
    M <- array(0, dim = c(nvars, nvars))
    if (i > 1) {
      if (is.null(groups)) {
        M <- outer(1:(nvars), 1:(nvars),
                   function(x, y) {
                     diff.x.y <- abs(beta.pen[y] - beta.pen[x])
                     diff.x.y[diff.x.y < tol] <- tol
                     ifelse(y > x, 1 / diff.x.y, 0)
                   })
      } else if (!all(is.na(groups))) {
        
        uniques <- unique(groups[!is.na(groups)])    
        
        for (t in 1:length(uniques)) {
          which.group <- which(groups == uniques[t])
          for (j in 1:(length(which.group) - 1)) {
            diff.x.y <- abs(beta.pen[which.group[j]] - beta.pen[which.group[j + 1]])
            
            if (diff.x.y < tol){diff.x.y <- tol}
            M[which.group[j], which.group[j+1]] <- 1 / diff.x.y
            M[which.group[j+1], which.group[j]] <- 1 / diff.x.y
          }
        }
      }
    }
    L <- lambda.fused * (diag(rowSums(M)) - M)
    beta <- ridge.fit(x, drop(y), drop(w), lambda = NULL, 
                      mult.factor = A + L, intercept = intercept)
    beta[abs(beta) < tol] <- tol
    
    if (all(abs(beta - prev) < tol)) {
      if (verbose) {cat("Converged at iteration: ", i, "\n")}
      break
    }
  }
  beta
}

fusedLogReg <- function(x, y, groups = NULL,
                        lambda.lasso = 0, lambda.fused = 0,
                        intercept = TRUE,
                        irls.maxiter = 30, irls.tol = 1e-10, 
                        fused.maxiter = 500, fused.tol = 1e-10) {
  
  y.working <- y
  
  nobs <- nrow(x)
  nvars <- ncol(x)
  w <- rep(0.5, nobs)
  len <- if (intercept) {nvars + 1} else {nvars}
  beta <- rep(1, len)
  for (i in 1:irls.maxiter) {
    prev <- beta
    
    init <- if (intercept) {prev[-1]} else {prev}
    
    beta.tmp <- fusedLassoRidge(x, y.working, w, groups = groups,
                                lambda.lasso = lambda.lasso, 
                                lambda.fused = lambda.fused,
                                maxiter = fused.maxiter,
                                intercept = FALSE,
                                tol = fused.tol, beta.init = init)
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
    
    if (all(abs(beta - prev) < fused.tol)) {
      cat("IRLS Converged at iteration: ", i, "\n")
      break
    }
  }
  beta
}

fusedMultinomLogReg <- function(x, y, groups = NULL,
                                lambda.lasso = 0, lambda.fused = 0,
                                intercept = TRUE,
                                irls.maxiter = 30, irls.tol = 1e-10, 
                                fused.maxiter = 500, fused.tol = 1e-10,
                                beta.init = NULL, groups.in = NULL) {
  
  y.f <- as.factor(y)
  classes <- levels(y.f)
  K <- length(classes)
  
  nobs <- nrow(x)
  nvars <- ncol(x)
  len <- if (intercept) {nvars + 1} else {nvars}
  if (!is.null(beta.init)) {
    stopifnot(all(dim(beta.init) == c(K, len)))
  }
  w <- rep(0.5, nobs)
  betas <- if(is.null(beta.init)) {array(1, dim = c(K, len))} else {beta.init}
  beta <- betas[1,]
  converged <- rep(FALSE, K)
  for (i in 1:irls.maxiter) {
    prev <- betas
    for (k in 1:K) {
      if (!converged[k]) { 
        y.working <- 1 * (y.f == classes[k])
        
        init <- if (intercept) {prev[k,-1]} else {prev[k,]}
        
        beta.tmp <- fusedLassoRidge(x, y.working, w, groups = groups,
                                    lambda.lasso = lambda.lasso, 
                                    lambda.fused = lambda.fused,
                                    maxiter = fused.maxiter,
                                    intercept = FALSE,
                                    tol = fused.tol, beta.init = init)
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
        
        if (all(abs(beta - prev[k,]) < fused.tol)) {
          converged[k] <- TRUE
        }
        
      }
    }
    if (all(converged)) {
      cat("IRLS Converged at iteration: ", i, "\n")
      break
    }
    cat("IRLS iter: ", i, "\n")
  }
  betas
}
