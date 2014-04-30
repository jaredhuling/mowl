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


groupFusedMultinomLogistic <- function(x, y, weights, groups = NULL,
                                           lambda = NULL, 
                                           nlambda = 21,
                                           intercept = TRUE,
                                           beta.init = NULL) {
  
  y.f <- as.factor(y)
  y <- as.numeric(levels(y.f)[y.f])
  classes <- levels(y.f)
  K <- length(classes)
  
  # EFLA options
  opts <- sllOpts()
  
  if (!is.null(lambda)) {
    stopifnot(length(lambda) == 3)
    lambda <- matrix(lambda, ncol = 3)
    colnames(lambda) <- c("lambda.lasso", "lambda.fused", "lambda.group")
    nlambda <- 1
  } else {
    # generate a lattice of 21 combinations
    # of tuning parameters based on a Unif Design
    lambda <- genLambdaLattice(nlambda, dim = 3) / ( sqrt(length(y)) )
    nlambda <- nrow(lambda)
  }

  
  
  nobs <- nrow(x)
  nvars <- ncol(x)
  len <- if (intercept) {nvars + 1} else {nvars}
  if (!is.null(beta.init)) {
    stopifnot(all(dim(beta.init) == c(K, len)))
  }
  betas <- if(is.null(beta.init)) {array(1, dim = c(K, len))} else {beta.init}
  beta <- betas[1,]
  beta.list <- funVal.list <- vector(mode = "list", length = nlambda)
  
  for (l in 1:nlambda) {
    current.lambdas <- lambda[l,]
    
    # solve the sparse group fused lasso multinomial logistic
    # regression problem 
    res <- fusedMultinomialLogistic(x, y, groups = groups, weights = weights,
                                    lambda = current.lambdas[1],
                                    lambda.fused = current.lambdas[2],
                                    lambda.group = current.lambdas[3],
                                    opts = opts)
    betas[,-1] <- res$beta
    betas[,1] <- res$intercept
    #betas <- res
    
    #initialize next one with current estimates
    opts$x0 <- t(res$beta)
    opts$c0 <- res$intercept
    
    attr(betas, "tuning.values") <- current.lambdas
    beta.list[[l]] <- betas
    funVal.list[[l]] <- res$funVal
    cat(l)
  } # end loop over tuning parameter combinations
  
  
  structure(list(coefficients = beta.list, lambda = lambda,
                 classes = classes, funVal = funVal.list), 
            class = "groupSparseFusedFit")
}




groupFusedMultinomLogReg2Stage <- function(x, y, weights, groups = NULL,
                                     lambda = NULL, nlambda.group = 50,
                                     nlambda.fused = 21,
                                     intercept = TRUE,
                                     irls.maxiter = 30, irls.tol = 1e-10,
                                     fused.maxiter = 500, fused.tol = 1e-10,
                                     beta.init = NULL, groups.in = NULL) {
  
  y <- as.factor(y)
  
  classes <- levels(y)
  
  # only do group lasso if group combinations not found
  if (is.null(groups.in)) {
    # select groups with Multinomial logistic regression 
    # with group lasso penalty
    gmlr <- groupMultinomLogReg(x, y, weights = weights, groups = groups, 
                                intercept = FALSE, nlambda = nlambda.group,
                                irls.maxiter = irls.maxiter, irls.tol = irls.tol)
    
    cat("group lasso stage complete \n")
    
    # find all unique combinations of selected groups
    # from group lasso penalty
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
  }
  
  # encourage sparsity with fused sparse lasso
  # with fused lasso penalty applied within groups
  fused.fit <- fusedMultinomLogRegSecond(x, y, weights = weights,
                                         groups = groups, intercept = intercept,
                                         nlambda = nlambda.fused, lambda = lambda,
                                         groups.in = groups.in)
  
  structure(list(coefficients = fused.fit$beta, lambda = fused.fit$lambda,
                 classes = classes, fused.iters = fused.fit$iters, groups.in = groups.in), 
              class = "groupSparseFusedFit2Stage")
}



fusedMultinomLogRegSecond <- function(x, y, weights = rep(1, nrow(x)), groups = NULL,
                                      nlambda = 21,
                                      lambda = NULL,
                                      intercept = TRUE,
                                      beta.init = NULL, groups.in = NULL) {
  
  y.f <- as.factor(y)
  y <- as.numeric(levels(y.f)[y.f])
  classes <- levels(y.f)
  K <- length(classes)
  G <- length(groups.in)
  
  # EFLA options
  opts <- sllOpts()
  
  if (!is.null(lambda)) {
    stopifnot(length(lambda) == 2)
    lambda <- matrix(lambda, ncol = 2)
    colnames(lambda) <- c("lambda.lasso", "lambda.fused")
    nlambda <- 1
  } else {
    # generate a lattice of 21 combinations
    # of tuning parameters based on a Unif Design
    lambda <- genLambdaLattice(nlambda) / sqrt(length(y))
  }
  
  nobs <- nrow(x)
  nvars <- ncol(x)
  len <- if (intercept) {nvars + 1} else {nvars}
  if (!is.null(beta.init)) {
    stopifnot(all(dim(beta.init) == c(K, len)))
  }
  w <- rep(0.5, nobs)
  betas <- if(is.null(beta.init)) {array(1, dim = c(K, len))} else {beta.init}
  beta <- betas[1,]
  beta.list <- iter.list <- vector(mode = "list", length = G)
  names(beta.list) <- as.character(1:G)
  grps <- sort(unique(groups))
  
  for (g in 1:G) {
    
    beta.tmp.list <- vector(mode = "list", length = nlambda)
    
    group.list <- nonzero.list <- vector(mode = "list", length = K)
     
    
    # set up lists (one element for each class)
    # to specify which groups to set to zero, as 
    # determined by group lasso
    for (k in 1:K) {
      which.groups <- grps[which(groups.in[[g]][k, ])]
      in.idx <- which(groups %in% which.groups | is.na(groups))
      group.list[[k]] <- groups[in.idx]
      nonzero.list[[k]] <- in.idx
    }

    
    for (l in 1:nlambda) {
      current.lambdas <- lambda[l,]
      
      # solve the sparse fused lasso multinomial logistic
      # regression problem with groups selected by group lasso 
      res <- fusedLassoMultinomLogisticStage2(x, y, 
                                              lambda.lasso = current.lambdas[1],
                                              lambda.fused = current.lambdas[2],
                                              group.list = group.list,
                                              nonzero.list = nonzero.list,
                                              opts = opts)
      #betas[,-1] <- res$beta
      #betas[,1] <- res$intercept
      betas <- res
      
      attr(betas, "tuning.values") <- current.lambdas
      beta.tmp.list[[l]] <- betas
      
    } # end loop over tuning parameter combinations
    cat("Group-lasso model", g,  "converged", "\n")
    beta.list[[g]] <- beta.tmp.list
    #iter.list[[g]] <- iter.tmp.list
  } # end loop over group-lasso values
  list(beta = beta.list, lambda = lambda) #, iters = iter.list)
}






