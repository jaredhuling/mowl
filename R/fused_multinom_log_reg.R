

fusedLassoMultinomLogisticStage2 <- function(x, y, group.list = NULL,
                                      nonzero.list = NULL,
                                      lambda.lasso = 0, lambda.fused = 0,
                                      intercept = TRUE, opts = NULL,
                                      irls.maxiter = 30, irls.tol = 1e-10, 
                                      beta.init = NULL, groups.in = NULL,
                                             verbose = FALSE) {
  
  y.f <- as.factor(y)
  classes <- levels(y.f)
  K <- length(classes)
  nobs <- nrow(x)
  nvars <- ncol(x)
  

  
  # if groups are given, get unique groups
  if (!is.null(group.list)) {
    for (k in 1:K) {
      groups <- group.list[[k]]
      unique.groups <- sort(unique(groups[!is.na(groups)]))
    }
  } else {
    group.list <- rep(list(rep(NA, nvars)), K)
  }
  
  nobs <- nrow(x)
  nvars <- ncol(x)
  len <- if (intercept) {nvars + 1} else {nvars}
  if (!is.null(beta.init)) {
    stopifnot(all(dim(beta.init) == c(K, len)))
  }
  
  betas <- if(is.null(beta.init)) {array(0, dim = c(K, len))} else {beta.init}
  beta <- betas[1,]
  z <- w <- vector(mode = "list", length = K)
  w[1:K] <- rep(list(rep(0.5, nobs)), K)
  converged <- rep(FALSE, K)
  for (i in 1:irls.maxiter) {
    prev <- betas
    for (k in 1:K) {
      if (!converged[k]) { 
        y.working <- 1 * (y.f == classes[k])
        
        if (i == 1) {
          z[[k]] <- 4 * (y.working - 0.5)
        }
        
        init <- if (intercept) {prev[k,-1]} else {prev[k,]}
        
        if (length(nonzero.list[[k]]) > 0) {
          print("z")
          print(z[[k]][1:10]); print("w")
          print(w[[k]][1:10])
          beta.tmp <- fusedlasso(x[,nonzero.list[[k]]], 
                                 z[[k]], w[[k]], 
                                 groups = group.list[[k]],
                                 lambda.lasso = lambda.lasso, 
                                 lambda.fused = lambda.fused, 
                                 family = "gaussian")
        
          if (intercept) {
            beta[nonzero.list[[k]]+1] <- beta.tmp$beta
          } else {
            beta[nonzero.list[[k]]] <- beta.tmp$beta
          }
        } else {
          beta <- rep(0, length(beta))
        }
            
        if (intercept) {
          xwb.tmp <- drop(x %*% beta[-1])
          #beta[1] <- mean( y.working - xwb.tmp)
          beta[1] <- beta.tmp$intercept
          xwb <- xwb.tmp + beta[1]
        } else {
          xwb <- drop(x %*% beta)
        }
        
        # update weights
        p <- 1 / (1 + exp(-xwb))
        w[[k]] <- p * (1 - p)
        w[[k]][which(w[[k]] < 1e-5)] <- 1e-5
        
        z[[k]] <- xwb + (y.working - p) / w[[k]]
        
        betas[k,] <- beta
        
        if (all(abs(beta - prev[k,]) < irls.tol)) {
          converged[k] <- TRUE
        }
        
      }
    }
    if (all(converged)) {
      cat("IRLS Converged at iteration: ", i, "\n")
      break
    }
    if (verbose) cat("IRLS iter: ", i, "\n")
  }
  betas
}



fusedMultinomialLogistic <- function(x, y, lambda, lambda.fused = 0,
                                     lambda.group = 0, 
                                     weights = rep(1, nrows(x)), 
                                     groups = NULL, 
                                     class.weights = NULL, opts=NULL) {
  
  sz <- dim(x)
  n <- sz[1]
  p <- sz[2]
  
  y.f <- as.factor(y)
  #if (is.factor(y)) {
  #  y <- levels(y)[y]
  #}
  classes <- levels(y.f)
  K <- length(classes)
  
  y.mat <- array(-1, dim = c(n, K))
  for (k in 1:K) {
    y.mat[y.f == classes[k],k] <- 1
  }
  
  betas <- array(NA, dim = c(K, p))
  intercepts <- numeric(K)
  
  stopifnot(lambda > 0)
  #if ( any(sort(unique(y)) != c(-1, 1)) ) {
  #  stop("y must be in {-1, 1}")
  #}
  opts.orig <- opts
  
  
  # if groups are given, get unique groups
  
  if (!is.null(groups)) {
    unique.groups <- vector(mode = "list", length = K)
    if (is.list(groups)) {
      if (length(groups) != K) {
        stop("Group list but have one element per class")
      }
      
      for (k in 1:K) {
        unique.groups[[k]] <- sort(unique(groups[!is.na(groups[[k]])]))
      }
    } else {
      unique.groups[1:K] <- rep(list(sort(unique(groups[!is.na(groups)]))), K)
      gr.list <- vector(mode = "list", length = K)
      gr.list[1:K] <- rep(list(groups), K)
      groups <- gr.list
    }
  }
  
  # run sllOpts to set default values (flags)
  opts <- sllOpts(opts.orig)
  if (lambda.group > 0) {
    opts$tol <- 1e-10
  }
  
  ## Set up options
  if (opts$nFlag != 0) {
    if (!is.null(opts$mu)) {
      mu <- opts$mu
      stopifnot(length(mu) == p)
    } else {
      mu <- colMeans(x)
    }
    
    if (opts$nFlag == 1) {
      if (!is.null(opts$nu)) {
        nu <- opts$nu
        stopifnot(length(nu) == p)
      } else {
        nu <- sqrt(colSums(x^2) / n)
      }
    }
    
    if (opts$nFlag == 2) {
      if (!is.null(opts$nu)) {
        nu <- opts$nu
        stopifnot(length(nu) == n)
      } else {
        nu <- sqrt(rowSums(x^2) / p)
      }
    }
    
    ## If some values of nu are small, it might
    ## be that the entries in a given row or col
    ## of x are all close to zero. For numerical 
    ## stability, we set these values to 1.
    ind.zero <- which(abs(nu) <= 1e-10)
    nu[ind.zero] <- 1
    
  }
  
  
  
  #code current y as 1 and -1
  #y.k <- 2 * (y.f == classes[k]) - 1
  
  ## Group & others
  
  # The parameter 'weight' contains the weight for each training sample.
  # See the definition of the problem above.
  # The summation of the weights for all the samples equals to 1.
  p.flag <- (y.mat == 1)
  if (!is.null(class.weights)) {
    sWeight <- class.weights
    
    if (length(class.weights) != 2 || class.weights[1] <= 0 || class.weights[2] <= 0) {
      stop("class weights must contain 2 positive values")
    }
    
    weight <- numeric(n)
    
    m1 <- sum(p.flag) * sWeight[1]
    m2 <- sum(!p.flag) * sWeight[2]
    weight[which(p.flag)] <- sWeight[1] / (m1 + m2)
    weight[which(!p.flag)] <- sWeight[2] / (m1 + m2)
  } else {
    weight <- array(1, dim = c(n, K)) / n
  }
  
  ## L2 norm regularization
  if (!is.null(opts$rsL2)) {
    rsL2 <- opts$rsL2
    if (rsL2 < 0) {
      stop("opts$rsL2 must be nonnegative")
    }
  } else {
    rsL2 <- 0
  }
  
  #m1 <- sum(weight[which(p.flag)])
  m1 <- colSums(p.flag * weight)
  m2 <- 1 - m1
  
  
  ## L1 norm regularization
  if (opts$rFlag != 0) {
    if (lambda < 0 || lambda > 1) {
      stop("opts.rFlag=1, and lambda should be in [0,1]")
    }
    
    ## we compute ATb for computing lambda_max, when the input lambda is a ratio
    b <- array(0, dim = c(n, K))
    b[which(p.flag)] <- m2
    b[which(!p.flag)] <- -m1
    b <- b * weight * weights
    #b <- b / n
    
    ## compute xTb
    if (opts$nFlag == 0) {
      xTb <- crossprod(x, b)
    } else if (opts$nFlag == 1) {
      xTb <- (crossprod(x, b) - colSums(b) * mu) / nu
    } else {
      invNu <- b / nu
      xTb <- crossprod(x, invNu) - colSums(invNu) * mu
    }
    
    lambda.max <- max(abs(xTb))
    lambda <- lambda * lambda.max
    
    if (is.null(opts$fusedPenalty)) {
      lambda2 <- 0
    } else {
      lambda2 <- lambda.max * lambda.fused
    }
    
    rsL2 <- rsL2 * lambda.max
    
  } else {
    lambda2 <- lambda.fused
  }
  
  ## initialize a starting point
  if (opts$init == 2) {
    beta <- array(0, dim = c(p, K))
    c <- log(m1 / m2)
  } else {
    if (!is.null(opts$x0)) {
      beta <- opts$x0
      if (any(dim(beta) != c(p, K))) {
        stop("initialization must be of dimension p x K")
      }
    } else {
      beta <- array(0, dim = c(p, K))
    }
    
    if (!is.null(opts$c0)) {
      c <- opts$c0
      if (length(c) != K) {
        stop("c0 must be of length K")
      }
    } else {
      c <- log(m1 / m2)
    }
  }
  
  ## Compute x beta
  if (opts$nFlag == 0) {
    xbeta <- x %*% beta
  } else if (opts$nFlag == 1) {
    invNu <- beta / nu
    mu.invNu <- as.double(crossprod(mu, invNu))
    xbeta <- x %*% invNu - rep(mu.invNu, n)
  } else {
    mubeta <- as.double(crossprod(mu, beta))
    xbeta <- (x %*% beta - rep(mubeta, n)) / nu
  }
  
  # restart the program for better efficiency
  # this is a newly added function
  if (is.null(opts$rStartNum)) {
    opts$rStartNum <- opts$maxIter
  } else {
    if (opts$rStartNum <= 0) {
      opts$rStartNum <- opts$maxIter
    }
  }
  
  ### THE MAIN CODE
  
  # z0 is the starting point for flsa
  z0 <- array(0, dim = c((p-1), K))
  
  ValueL <- funVal <- numeric(opts$maxIter)
  
  ## The Armijo Goldstein line search scheme + accelerated gradient descent
  if (opts$mFlag == 0 && opts$lFlag == 0) {
    
    # this flag tests whether the gradient step only changes a little
    bFlag <- 0 
    
    # the intial guess of the Lipschitz continuous gradient
    L <- 1 / n + rsL2
    
    # the product between weight and y
    #weighty <- weight * y.k
    
    #initial assignments
    betap <- beta
    xbetap <- xbeta
    betabetap <- numeric(p)
    cp <- c; ccp <- numeric(K)
    
    alphap <- 0; alpha <- 1
    
    for (iterStep in 1:opts$maxIter) {
      
      diff.prev <- -n
      
      ## --------------------------- step 1 ---------------------------
      ##      compute search point s based on xp and x (with beta)
      bet <- (alphap - 1) / alpha
      s <- beta + bet * betabetap
      sc <- c + bet * ccp
      
      ## --------------------------- step 2 ---------------------------
      ##  line search for L and compute the new approximate solution x
      
      # compute xs = x * s
      xs <- xbeta + bet * (xbeta - xbetap)
      
      #aa <- -y.k * (xs + sc)
      aa <- (xs + rep(sc, each = n))
      
      # fun.s is the logistic loss at the search point
      bb <- pmax(- y.mat * aa, 0)
      #fun.s <- as.double( crossprod(weight, (log(exp(-bb) + exp(aa - bb)) + bb)) ) + 
      #  ( rsL2 / 2 ) * as.double(crossprod(s))
      
      
      # compute prob=[p_1;p_2;...;p_n]
      #prob <- 1 / (1 + exp(aa))
      prob <- exp(aa)
      
      
      rSp <- rowSums(prob)
      
      fun.s <- -sum(weights * (rowSums(((y.mat + 1) / 2) * aa) - log( rSp ))) / n + 
        ( rsL2 / 2 ) * sum(as.double(crossprod(s)))
      
      prob <- prob / rSp
      
      #b <- -weighty * (1 - prob)
      b <- -((y.mat+1)/2 - prob) * weight * weights
      
      #the gradient of c
      #gc <- (colSums(y.mat+1)/(2) - 1) / n
      gc <- colSums(b)
      
      #  should be sum i=1:n { sum k=1:K {y_i^(k)} - p_ij} 
      
      #compute g= xT b, the gradient of beta
      if (opts$nFlag == 0) {
        g <- crossprod(x, b)
      } else if (opts$nFlag == 1) {
        g <- (crossprod(x, b) - colSums(b) * mu) / nu
      } else {
        invNu <- b / nu
        g <- crossprod(x, invNu) - colSums(invNu) * mu
      }
      #add the squared L2 norm regularization term
      g <- g + rsL2 * s
      
      #assignments
      betap <- beta
      xbetap <- xbeta
      cp <- c
      
      while (TRUE) {
        # let s walk in a step in the antigradient of s to get v
        # and then do the Lq/L1-norm regularized projection
        v <- s - g / L
        c <- sc - gc / L
        
        for (k in 1:K) {
          if (is.null(groups)) {
            res <- flsa(v[, k], z0[, k], lambda / L, lambda2 / L, p,
                        1000, 1e-9, 1, 6)
            beta[, k] <- res[[1]]
            z0[, k] <- res[[2]]
            infor <- res[[3]]
          } else {
            
            if (any(is.na(groups[[k]]))) {
              ## don't apply fused lasso penalty
              ## to variables with group == NA 
              gr.idx <- which(is.na(groups[[k]]))
              gr.p <- length(gr.idx)
              if (any(gr.idx == 1)) {
                gr.idx.z <- gr.idx[gr.idx != 1] - 1
              } else {
                gr.idx.z <- gr.idx[-gr.p]
              }
              
              res <- flsa(v[gr.idx, k], z0[gr.idx.z, k], lambda / L, 0, gr.p,
                          1000, 1e-9, 1, 6)            
              
              beta[gr.idx, k] <- res[[1]]
              z0[gr.idx.z, k] <- res[[2]]
              infor <- res[[3]]
            }
            
            for (t in 1:length(unique.groups[[k]])) {
              gr.idx <- which(groups[[k]] == unique.groups[[k]][t])
              gr.p <- length(gr.idx)
              if (any(gr.idx == 1)) {
                gr.idx.z <- gr.idx[gr.idx != 1] - 1
              } else {
                gr.idx.z <- gr.idx[-gr.p]
              }
              
              res <- flsa(v[gr.idx, k], z0[gr.idx.z, k], lambda / L, lambda2 / L, gr.p,
                          1000, 1e-9, 1, 6)
              
              if (lambda.group > 0) {
                ## 2nd Projection:
                ## argmin_w { 0.5 \|w - w_1\|_2^2
                ##          + lambda_3 * \|w_1\|_2 }
                ## This is a simple thresholding:
                ##    w_2 = max(\|w_1\|_2 - \lambda_3, 0)/\|w_1\|_2 * w_1
                nm = norm(res[[1]], type = "2")
                if (nm == 0) {
                  newbeta = numeric(length(res[[1]]))
                } else {
                  #apply soft thresholding, adjust penalty for size of group
                  newbeta = pmax(nm - lambda.group * sqrt(gr.p), 0) / nm * res[[1]]
                }
                end
              } else {
                newbeta <- res[[1]]
              }
              
              beta[gr.idx, k] <- newbeta
              z0[gr.idx.z, k] <- res[[2]]
              infor <- res[[3]]
            }
          }
        } #end loop over classes
        
        # the difference between the new approximate 
        # solution x and the search point s
        v <- beta - s
        
        ## Compute x beta
        if (opts$nFlag == 0) {
          xbeta <- x %*% beta
        } else if (opts$nFlag == 1) {
          invNu <- beta / nu
          mu.invNu <- as.double(crossprod(mu, invNu))
          xbeta <- x %*% invNu - rep(mu.invNu, n)
        } else {
          mubeta <- as.double(crossprod(mu, beta))
          xbeta <- (x %*% beta - rep(mubeta, n)) / nu
        }
        
        #aa <- -y.k * (xbeta + c)
        aa <- (xbeta + rep(c, each = n))
        
        # fun.beta is the logistic loss at the new approximate solution
        bb <- pmax(- y.mat * aa, 0)
        
        
        fun.beta <- -sum(weights * (rowSums(((y.mat + 1) / 2) * aa) - log( rowSums( exp(aa) )) ) ) / n + 
          ( rsL2 / 2 ) * sum(as.double(crossprod(beta)))
        #in case of super bad likelihood, just make it reasonably big
        if (is.nan(fun.beta) | fun.beta == Inf) {fun.beta <- 1e10}  
        
        r.sum <- norm(v, type = "F") ^ 2 / 2 + sum((c - sc)^2) / 2
        fzp.gamma <- fun.s + sum(sum(v * g)) + L * r.sum + sum((c - sc) * gc)

        if (r.sum <= 1e-18) {
          #this shows that the gradient step makes little improvement
          bFlag <- 1
          break
        }
        
        ## the condition is fun.beta <= fun.s + v'* g + c * gc
        ##                           + L/2 * (v'*v + (c-sc)^2 )
        
        if (fun.beta <= fzp.gamma ) { #| funval.s - funval < diff.prev
          break
        } else {
          #L <- max(2 * L, (fun.beta * L) / fzp.gamma)
          L <- 2 * L
        }
        
        #diff.prev <- funval.s - funval
        
      } # end while loop
      
      ## --------------------------- step 3 ---------------------------
      ##      update alpha and alphap, and check for convergence
      alphap <- alpha
      alpha <- (1 + sqrt(4 * alpha * alpha + 1)) / 2
      
      ## store values for L
      ValueL[iterStep] <- L
      
      betabetap <- beta - betap
      ccp <- c - cp
      
      
      # evaluate fused and group lasso 
      # penalty-terms
      if (!is.null(groups)) {
        fused.pen <- group.pen <- 0
        for (k in 1:K) {
          for (t in 1:length(unique.groups[[k]])) {
            gr.idx <- which(groups[[k]] == unique.groups[[k]][t])
            gr.p <- length(gr.idx)
            if (gr.p > 1) {
              fused.pen <- fused.pen + sum(abs(beta[gr.idx[2:(gr.p)], k] - beta[gr.idx[1:(gr.p - 1)], k]))
              group.pen <- group.pen + sqrt(sum(beta[gr.idx, k] ^ 2) * gr.p)
            }
          }
        }
        pens <- lambda2 * fused.pen + lambda.group * group.pen
      } else {
        pens <- lambda2 * sum(abs(beta[2:p] - beta[1:(p-1)]))
      }
      
      funVal[iterStep] <- fun.beta + lambda * sum(abs(beta)) + pens
      
      if (bFlag) {
        break
      }
      
      tf <- opts$tFlag
      
      if (tf == 0) {
        if (iterStep > 2) {
          if (abs( funVal[iterStep] - funVal[iterStep - 1] ) <= opts$tol) {
            break
          }
        }
      } else if (tf == 1) {
        if (iterStep > 2) {
          if (abs( funVal[iterStep] - funVal[iterStep - 1] ) <= 
                opts$tol * funVal[iterStep - 1]) {
            break
          }
        }
      } else if (tf == 2) {
        if (funVal[iterStep] <= opts$tol) {
          break
        }
      } else if (tf == 3) {
        norm.bbp <- sqrt(as.double(crossprod(bbp)))
        if (norm.bbp <= opts$tol) {
          break
        }
      } else if (tf == 4) {
        norm.bp <- sqrt(as.double(crossprod(bp)))
        norm.bbp <- sqrt(as.double(crossprod(bbp)))
        if (norm.bbp <= opts$tol * max(norm.bp)) {
          break
        }
      } else if (tf == 5) {
        if (iterStep >= opts$maxIter) {
          break
        }
      }
      
      # restart the program every opts$rStartNum
      if ((iterStep %% opts$rStartNum == 0)) {
        alphap <- 0; alpha <- 1
        betap <- beta; xbetap <- xbeta; L <- 1
        betabetap <- numeric(p)
      }
      
    } #end of iterStep loop
    
    funVal <- funVal[1:iterStep]
    ValueL <- ValueL[1:iterStep]
    
  } else {
    stop("This function only supports opts.mFlag=0 & opts.lFlag=0")
  }
  
  
  list(beta = t(beta), intercept = c, funVal = funVal, ValueL = ValueL, bFlag = bFlag)
}






fusedLassoMultinomLogisticStage2_false <- function(x, y, lambda, lambda.fused = NULL,
                                                   group.list = NULL, 
                                                   nonzero.list = NULL,
                                                   class.weights = NULL, opts=NULL) {
  
  sz <- dim(x)
  n <- sz[1]
  p <- sz[2]
  
  y.f <- as.factor(y)
  #if (is.factor(y)) {
  #  y <- levels(y)[y]
  #}
  classes <- levels(y.f)
  K <- length(classes)
  
  betas <- array(0, dim = c(K, p))
  intercepts <- numeric(K)
  
  stopifnot(lambda > 0)
  #if ( any(sort(unique(y)) != c(-1, 1)) ) {
  #  stop("y must be in {-1, 1}")
  #}
  opts.orig <- opts
  
  for (k in 1:K) {
    
    # if groups are given, get unique groups
    if (!is.null(group.list)) {
      groups <- group.list[[k]]
      unique.groups <- sort(unique(groups[!is.na(groups)]))
      x.g <- x[,nonzero.list[[k]]]
      sz <- dim(x.g)
      n <- sz[1]
      p <- sz[2]
    } else {
      x.g <- x
    }
    
    # run sllOpts to set default values (flags)
    opts <- sllOpts(opts.orig)
    
    ## Set up options
    if (opts$nFlag != 0) {
      if (!is.null(opts$mu)) {
        mu <- opts$mu
        stopifnot(length(mu) == p)
      } else {
        mu <- colMeans(x.g)
      }
      
      if (opts$nFlag == 1) {
        if (!is.null(opts$nu)) {
          nu <- opts$nu
          stopifnot(length(nu) == p)
        } else {
          nu <- sqrt(colSums(x.g^2) / n)
        }
      }
      
      if (opts$nFlag == 2) {
        if (!is.null(opts$nu)) {
          nu <- opts$nu
          stopifnot(length(nu) == n)
        } else {
          nu <- sqrt(rowSums(x.g^2) / p)
        }
      }
      
      ## If some values of nu are small, it might
      ## be that the entries in a given row or col
      ## of x are all close to zero. For numerical 
      ## stability, we set these values to 1.
      ind.zero <- which(abs(nu) <= 1e-10)
      nu[ind.zero] <- 1
      
    }
    
    
    
    #code current y as 1 and -1
    y.k <- 2 * (y.f == classes[k]) - 1
    
    ## Group & others
    
    # The parameter 'weight' contains the weight for each training sample.
    # See the definition of the problem above.
    # The summation of the weights for all the samples equals to 1.
    p.flag <- (y.k == 1)
    if (!is.null(class.weights)) {
      sWeight <- class.weights
      
      if (length(class.weights) != 2 || class.weights[1] <= 0 || class.weights[2] <= 0) {
        stop("class weights must contain 2 positive values")
      }
      
      weight <- numeric(n)
      
      m1 <- sum(p.flag) * sWeight[1]
      m2 <- sum(!p.flag) * sWeight[2]
      weight[which(p.flag)] <- sWeight[1] / (m1 + m2)
      weight[which(!p.flag)] <- sWeight[2] / (m1 + m2)
    } else {
      weight <- rep(1, n) / n
    }
    
    ## L2 norm regularization
    if (!is.null(opts$rsL2)) {
      rsL2 <- opts$rsL2
      if (rsL2 < 0) {
        stop("opts$rsL2 must be nonnegative")
      }
    } else {
      rsL2 <- 0
    }
    
    m1 <- sum(weight[which(p.flag)])
    m2 <- 1 - m1
    
    ## L1 norm regularization
    if (opts$rFlag != 0) {
      if (lambda < 0 || lambda > 1) {
        stop("opts.rFlag=1, and lambda should be in [0,1]")
      }
      
      ## we compute ATb for computing lambda_max, when the input lambda is a ratio
      b <- numeric(n)
      b[which(p.flag)] <- m2
      b[which(!p.flag)] <- -m1
      b <- b * weight
      
      ## compute xTb
      if (opts$nFlag == 0) {
        xTb <- crossprod(x.g, b)
      } else if (opts$nFlag == 1) {
        xTb <- (crossprod(x.g, b) - sum(b) * mu) / nu
      } else {
        invNu <- b / nu
        xTb <- crossprod(x.g, invNu) - sum(invNu) * mu
      }
      
      lambda.max <- max(abs(xTb))
      lambda <- lambda * lambda.max
      
      if (is.null(lambda.fused)) {
        lambda2 <- 0
      } else {
        lambda2 <- lambda.max * lambda.fused
      }
      
      rsL2 <- rsL2 * lambda.max
      
    } else {
      if (is.null(lambda.fused)) {
        lambda2 <- 0
      } else {
        lambda2 <- lambda.fused
      }
    }
    
    ## initialize a starting point
    if (opts$init == 2) {
      beta <- numeric(p)
      c <- log(m1 / m2)
    } else {
      if (!is.null(opts$x0)) {
        beta <- opts$x0
        if (length(beta) != p) {
          stop("initialization must be of length p")
        }
      } else {
        beta <- numeric(p)
      }
      
      if (!is.null(opts$c0)) {
        c <- opts$c0
      } else {
        c <- log(m1 / m2)
      }
    }
    
    ## Compute x beta
    if (opts$nFlag == 0) {
      xbeta <- x.g %*% beta
    } else if (opts$nFlag == 1) {
      invNu <- beta / nu
      mu.invNu <- as.double(crossprod(mu, invNu))
      xbeta <- x.g %*% invNu - rep(mu.invNu, n)
    } else {
      mubeta <- as.double(crossprod(mu, beta))
      xbeta <- (x.g %*% beta - rep(mubeta, n)) / nu
    }
    
    # restart the program for better efficiency
    # this is a newly added function
    if (is.null(opts$rStartNum)) {
      opts$rStartNum <- opts$maxIter
    } else {
      if (opts$rStartNum <= 0) {
        opts$rStartNum <- opts$maxIter
      }
    }
    
    ### THE MAIN CODE
    
    # z0 is the starting point for flsa
    z0 <- numeric(p-1)
    
    ValueL <- funVal <- numeric(opts$maxIter)
    
    ## The Armijo Goldstein line search scheme + accelerated gradient descent
    if (opts$mFlag == 0 && opts$lFlag == 0) {
      
      # this flag tests whether the gradient step only changes a little
      bFlag <- 0 
      
      # the intial guess of the Lipschitz continuous gradient
      L <- 1 / n + rsL2
      
      # the product between weight and y
      weighty <- weight * y.k
      
      #initial assignments
      betap <- beta
      xbetap <- xbeta
      betabetap <- numeric(p)
      cp <- c; ccp <- 0
      
      alphap <- 0; alpha <- 1
      
      for (iterStep in 1:opts$maxIter) {
        
        ## --------------------------- step 1 ---------------------------
        ##      compute search point s based on xp and x (with beta)
        bet <- (alphap - 1) / alpha
        s <- beta + bet * betabetap
        sc <- c + bet * ccp
        
        ## --------------------------- step 2 ---------------------------
        ##  line search for L and compute the new approximate solution x
        
        # compute xs = x * s
        xs <- xbeta + bet * (xbeta - xbetap)
        
        aa <- -y.k * (xs + sc)
        
        # fun.s is the logistic loss at the search point
        bb <- pmax(aa, 0)
        fun.s <- as.double( crossprod(weight, (log(exp(-bb) + exp(aa - bb)) + bb)) ) + 
          ( rsL2 / 2 ) * as.double(crossprod(s))
        
        # compute prob=[p_1;p_2;...;p_n]
        prob <- 1 / (1 + exp(aa))
        
        b <- -weighty * (1 - prob)
        
        #the gradient of c
        gc <- sum(b) 
        
        #compute g= xT b, the gradient of beta
        if (opts$nFlag == 0) {
          g <- crossprod(x.g, b)
        } else if (opts$nFlag == 1) {
          g <- (crossprod(x.g, b) - sum(b) * mu) / nu
        } else {
          invNu <- b / nu
          g <- crossprod(x.g, invNu) - sum(invNu) * mu
        }
        #add the squared L2 norm regularization term
        g <- g + rsL2 * s
        
        #assignments
        betap <- beta
        xbetap <- xbeta
        cp <- c
        
        while (TRUE) {
          # let s walk in a step in the antigradient of s to get v
          # and then do the Lq/L1-norm regularized projection
          v <- s - g / L
          c <- sc - gc / L
          
          if (is.null(groups)) {
            res <- flsa(v, z0, lambda / L, lambda2 / L, p,
                        1000, 1e-8, 1, 6)
            beta <- res[[1]]
            z0 <- res[[2]]
            infor <- res[[3]]
          } else {
            
            if (any(is.na(groups))) {
              ## don't apply fused lasso penalty
              ## to variables with group == NA 
              gr.idx <- which(is.na(groups))
              gr.p <- length(gr.idx)
              if (gr.p < 1) {
                next
              }
              if (any(gr.idx == 1)) {
                gr.idx.z <- gr.idx[gr.idx != 1] - 1
              } else {
                gr.idx.z <- gr.idx[-gr.p]
              }
              
              res <- flsa(v[gr.idx], z0[gr.idx.z], lambda / L, 0, gr.p,
                          1000, 1e-8, 1, 6)
              beta[gr.idx] <- res[[1]]
              z0[gr.idx.z] <- res[[2]]
              infor <- res[[3]]
            }
            
            for (t in 1:length(unique.groups)) {
              gr.idx <- which(groups == unique.groups[t])
              gr.p <- length(gr.idx)
              if (gr.p < 1) {
                next
              }
              if (any(gr.idx == 1)) {
                gr.idx.z <- gr.idx[gr.idx != 1] - 1
              } else {
                gr.idx.z <- gr.idx[-gr.p]
              }
              
              res <- flsa(v[gr.idx], z0[gr.idx.z], lambda / L, lambda2 / L, gr.p,
                          1000, 1e-8, 1, 6)
              beta[gr.idx] <- res[[1]]
              z0[gr.idx.z] <- res[[2]]
              infor <- res[[3]]
            }
          }
          
          # the difference between the new approximate 
          # solution x and the search point s
          v <- beta - s
          
          ## Compute x beta
          if (opts$nFlag == 0) {
            xbeta <- x.g %*% beta
          } else if (opts$nFlag == 1) {
            invNu <- beta / nu
            mu.invNu <- as.double(crossprod(mu, invNu))
            xbeta <- x.g %*% invNu - rep(mu.invNu, n)
          } else {
            mubeta <- as.double(crossprod(mu, beta))
            xbeta <- (x.g %*% beta - rep(mubeta, n)) / nu
          }
          
          aa <- -y.k * (xbeta + c)
          
          # fun.beta is the logistic loss at the new approximate solution
          bb <- pmax(aa, 0)
          
          fun.beta <- as.double( crossprod(weight, (log(exp(-bb) + exp(aa - bb)) + bb)) ) + 
            ( rsL2 / 2 ) * as.double(crossprod(beta))
          
          r.sum <- (as.double(crossprod(v)) + (c - sc)^2) / 2
          l.sum <- fun.beta - fun.s - as.double(crossprod(v, g)) - (c - sc) * gc
          
          
          if (r.sum <= 1e-18) {
            #this shows that the gradient step makes little improvement
            bFlag <- 1
            break
          }
          
          ## the condition is fun.beta <= fun.s + v'* g + c * gc
          ##                           + L/2 * (v'*v + (c-sc)^2 )
          
          if (l.sum <= r.sum * L) {
            break
          } else {
            L <- max(2 * L, l.sum / r.sum)
          }
          
        } # end while loop
        
        ## --------------------------- step 3 ---------------------------
        ##      update alpha and alphap, and check for convergence
        alphap <- alpha
        alpha <- (1 + sqrt(4 * alpha * alpha + 1)) / 2
        
        ## store values for L
        ValueL[iterStep] <- L
        
        betabetap <- beta - betap
        ccp <- c - cp
        
        funVal[iterStep] <- fun.beta + lambda * sum(abs(beta)) +
          lambda2 * sum(abs(beta[2:p] - beta[1:(p-1)]))
        
        if (bFlag) {
          break
        }
        
        tf <- opts$tFlag
        
        if (tf == 0) {
          if (iterStep > 2) {
            if (abs( funVal[iterStep] - funVal[iterStep - 1] ) <= opts$tol) {
              break
            }
          }
        } else if (tf == 1) {
          if (iterStep > 2) {
            if (abs( funVal[iterStep] - funVal[iterStep - 1] ) <= 
                  opts$tol * funVal[iterStep - 1]) {
              break
            }
          }
        } else if (tf == 2) {
          if (funVal[iterStep] <= opts$tol) {
            break
          }
        } else if (tf == 3) {
          norm.bbp <- sqrt(as.double(crossprod(bbp)))
          if (norm.bbp <= opts$tol) {
            break
          }
        } else if (tf == 4) {
          norm.bp <- sqrt(as.double(crossprod(bp)))
          norm.bbp <- sqrt(as.double(crossprod(bbp)))
          if (norm.bbp <= opts$tol * max(norm.bp)) {
            break
          }
        } else if (tf == 5) {
          if (iterStep >= opts$maxIter) {
            break
          }
        }
        
        # restart the program every opts$rStartNum
        if ((iterStep %% opts$rStartNum == 0)) {
          alphap <- 0; alpha <- 1
          betap <- beta; xbetap <- xbeta; L <- 1
          betabetap <- numeric(p)
        }
        
      } #end of iterStep loop
      
      funVal <- funVal[1:iterStep]
      ValueL <- ValueL[1:iterStep]
      
    } else {
      stop("This function only supports opts.mFlag=0 & opts.lFlag=0")
    }
    betas[k, nonzero.list[[k]]] <- drop(beta)
    intercepts[k] <- c
  }
  
  list(beta = betas, intercept = intercepts, funVal = funVal, ValueL = ValueL)
}