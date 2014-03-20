simulateOwlData <- function(n, p, rules, true.beta) {
  
  x <- matrix(rnorm(n * p), n, p)
  #x <- matrix(rbinom(n * p, 1, 0.05), n, p)
  
  interaction <- if(length(rules[[1]]) == 4) {TRUE} else {FALSE}
  
  
  if (interaction) {
    x <- genInteractions(x)
  }
  
  A <- as.factor(sample(3, n, replace = TRUE))
  
  treatment.effects <- genTreatmentEffects(x, A, treatment.rules)
  
  y <- rnorm(n) + rowSums(treatment.effects) + x[,3:(3 + length(true.beta) - 1)] %*% true.beta
  y <- y + abs(min(y))
  
  oracle <- genOracleTreatments(x, treatment.rules)
  ret <- list(x = x, y = y, A = A, oracle = oracle, treatment.effects = treatment.effects)
  class(ret) <- "sim.owl.data"
  ret
}

genInteractions <- function(x) {
  nobs <- nrow(x); nvars <- ncol(x)
  x.int <- array(0, dim = c(nobs, choose(nvars, 2)))
  k <- 0
  for (i in 1:(nvars - 1)) {
    for (j in (i + 1):nvars) {
      k <- k + 1
      x.int[,k] <- x[,i] * x[,j]
    }
  }
  cbind(x, x.int)
}