simulateOwlData <- function(n, p, rules, true.beta, interaction, sd.x = 1, sd.y = 1,
                            binary.outcome = FALSE, num.factors = 0L, factor.levels = rep(3, num.factors)) {
  
  num.contin <- p - num.factors
  num.factors <- as.integer(num.factors)
  n <- as.integer(n)
  p <- as.integer(p)
  x <- matrix(rnorm(n * num.contin, sd = sd.x), n, num.contin)
  colnames(x) <- paste("X", 1:p, sep = "")
  #x <- matrix(rbinom(n * p, 1, 0.05), n, p)
  x.factors <- sapply(factor.levels, function(x) sample.int(x, n, replace = TRUE))
  x <- cbind(x, x.factors)
    
  if (interaction) {
    x <- model.matrix(~ -1 + . + .*., data = x)
  } else {
    x <- model.matrix(~ -1 + ., data = x)
  }
  
  A <- as.factor(sample(3, n, replace = TRUE))
  
  treatment.effects <- genTreatmentEffects(x, A, rules)
  
  if (binary.outcome) {
    log.p.ratio <- rowSums(treatment.effects) + x[,3:(3 + length(true.beta) - 1)] %*% true.beta
    prob.y.1 <- 1 / (1 + exp(-log.p.ratio))
    y <- rbinom(n, 1, prob = prob.y.1)
  } else {
    y <- rnorm(n, sd = sd.y) + rowSums(treatment.effects) + x[,3:(3 + length(true.beta) - 1)] %*% true.beta
    y <- y + abs(min(y))
  }
  
  oracle <- genOracleTreatments(x, rules)
  ret <- list(x = x, y = y, A = A, oracle = oracle, treatment.effects = treatment.effects)
  class(ret) <- "sim.owl.data"
  ret
}

simulateOwlData2 <- function(n, p, rules, true.beta, interaction, sd.x = 1, sd.y = 1,
                            binary.outcome = FALSE, num.factors = 0L, factor.levels = rep(3, num.factors)) {
  
  num.contin <- p - num.factors
  num.factors <- as.integer(num.factors)
  n <- as.integer(n)
  p <- as.integer(p)
  x <- matrix(rnorm(n * num.contin, sd = sd.x), n, num.contin)
  
  #x <- matrix(rbinom(n * p, 1, 0.05), n, p)
  x.factors <- sapply(factor.levels, function(x) as.factor(sample.int(x, n, replace = TRUE)))
  x <- data.frame(x, x.factors)
  colnames(x) <- paste("X", 1:p, sep = "")
  
  if (interaction) {
    x <- model.matrix(~ -1 + . + .*., data = x)
  } else {
    x <- model.matrix(~ -1 + ., data = x)
  }
  
  A <- as.factor(sample(3, n, replace = TRUE))
  
  treatment.effects <- genTreatmentEffects2(x, A, rules)
  
  if (binary.outcome) {
    log.p.ratio <- rowSums(treatment.effects) + x[,3:(3 + length(true.beta) - 1)] %*% true.beta
    prob.y.1 <- 1 / (1 + exp(-log.p.ratio))
    y <- rbinom(n, 1, prob = prob.y.1)
  } else {
    y <- rnorm(n, sd = sd.y) + rowSums(treatment.effects) + x[,3:(3 + length(true.beta) - 1)] %*% true.beta
    y <- y + abs(min(y))
  }
  
  oracle <- genOracleTreatments2(x, rules)
  ret <- list(x = x, y = y, A = A, oracle = oracle, treatment.effects = treatment.effects)
  class(ret) <- "sim.owl.data"
  ret
}


genInteractions <- function(x) {
  nobs <- nrow(x); nvars <- ncol(x)
  x.int <- array(0, dim = c(nobs, choose(nvars, 2)))
  colnames(x.int) <- 1:choose(nvars, 2)
  k <- 0
  for (i in 1:(nvars - 1)) {
    for (j in (i + 1):nvars) {
      k <- k + 1
      x.int[,k] <- x[,i] * x[,j]
      colnames(x.int)[k] <- paste(colnames(x)[i], colnames(x)[j], sep = ":")
    }
  }
  cbind(x, x.int)
}