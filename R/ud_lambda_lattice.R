
genLambdaLattice <- function(nlambda = 21) {
  
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
  ud <- matrix(NA, nrow = n.lambda, ncol = 2)
  for(s in 1:2){
    ud[,s] <- ((1:n.lambda) * h[s] - 0.5) / n.lambda  
  }
  ud <- ud - floor(ud)
  colnames(ud) <- c("lambda.lasso", "lambda.fused")
  ud
}