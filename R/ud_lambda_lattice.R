
genLambdaLattice <- function(nlambda = 21, dim = 2) {
  
  if (dim == 2) {
    # set up possible uniform designs
    lbd <- c(8, 13, 21, 34, 55, 89, 144, 233, 377, 610)
    h2 <- c(5, 8, 13, 21, 34, 55, 89, 144, 233, 377)
    h <- data.frame(rbind(lbd, h2))
    rownames(h) <- c("n", "h2")
    
    
    if(!any(nlambda==lbd)){
      nlambda <- 21
      cat("Invalid nlambda, defaulting to ", nlambda, 
          " points of lambdas are sampled by uniform design.\n")
    }
    
    
    h=h[-1,which(h[1,] == nlambda)]
    h=c(1,h)
    
    # set up lattice
    ud <- matrix(NA, nrow = nlambda, ncol = 2)
    for(s in 1:2){
      ud[,s] <- ((1:nlambda) * h[s] - 0.5) / nlambda  
    }
    ud <- ud - floor(ud)
    colnames(ud) <- c("lambda.lasso", "lambda.fused")
  } else if (dim == 3) {
    
    lbd <- c(35, 101, 135, 185, 266, 418, 579, 828, 1010)
    h2 <- c(11,  40,  29,  26,  27,  90,  63, 285,  140)
    h3 <- c(16,  85,  42,  64,  69, 130, 169, 358,  237)
    h <- data.frame(rbind(lbd, h2, h3))
    rownames(h) <- c("n", "h2", "h3")
    
    
    if(!any(nlambda==lbd)){
      nlambda <- 35
      cat("Invalid nlambda, defaulting to ", nlambda, 
          " points of lambdas are sampled by uniform design.\n nlambda must be in: ", lbd)
    }
    
    
    h=h[-1,which(h[1,] == nlambda)]
    h=c(1,h)
    
    # set up lattice
    ud <- matrix(NA, nrow = nlambda, ncol = 3)
    for(s in 1:3){
      ud[,s] <- ((1:nlambda) * h[s] - 0.5) / nlambda  
    }
    ud <- ud - floor(ud)
    colnames(ud) <- c("lambda.lasso", "lambda.fused", "lambda.group")
    
  } else {
    stop("only dim = 2 and 3 supported")
  }
  ud
}