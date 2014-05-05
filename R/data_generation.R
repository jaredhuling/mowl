simulateOwlData <- function(n, p, rules, true.beta, interaction, sd.x = 1, sd.y = 1,
                            outcome.type = c("numeric", "binary", "misspecified.binary"), 
                            num.factors = 0L, factor.levels = rep(3, num.factors),
                            misspecified = FALSE) {
  
  outcome.type <- match.arg(outcome.type)
  num.contin <- p - num.factors
  num.factors <- as.integer(num.factors)
  n <- as.integer(n)
  p <- as.integer(p)
  x <- matrix(rnorm(n * num.contin, sd = sd.x), n, num.contin)
  colnames(x) <- paste("X", 1:p, sep = "")
  #x <- matrix(rbinom(n * p, 1, 0.05), n, p)
  if (num.factors > 0){
    x.factors <- sapply(factor.levels, function(x) sample.int(x, n, replace = TRUE))
    x <- data.frame(x, x.factors)
  } else {
    x <- data.frame(x)
  }
  if (interaction) {
    x <- model.matrix(~ -1 + . + .*., data = x)
  } else {
    x <- model.matrix(~ -1 + ., data = x)
  }
  
  A <- as.factor(sample(3, n, replace = TRUE))
  
  treatment.effects <- genTreatmentEffects(x, A, rules)
  
  if (misspecified) {
    y.main.eff <- rnorm(n * length(true.beta), ncol = length(true.beta)) %*% true.beta
  } else {
    y.main.eff <- x[,3:(3 + length(true.beta) - 1)] %*% true.beta
  }
  
  
  if (outcome.type == "binary") {
    log.p.ratio <- rowSums(treatment.effects) + y.main.eff
    prob.y.1 <- 1 / (1 + exp(-log.p.ratio))
    y <- rbinom(n, 1, prob = prob.y.1)
  } else if (outcome.type == "numeric") {
    y <- rnorm(n, sd = sd.y) + rowSums(treatment.effects) + y.main.eff
    y <- y + abs(min(y))
  } else if (outcome.type == "misspecified.binary"){
    exp.p.ratio <- rowSums(treatment.effects) + y.main.eff + 1
    prob.y.1 <- log(exp.p.ratio + sign(min(exp.p.ratio)) * min(exp.p.ratio) + 1) / (1 + log(exp.p.ratio + sign(min(exp.p.ratio)) * min(exp.p.ratio) + 1))
    y <- rbinom(n, 1, prob = prob.y.1)
  }
  
  oracle <- genOracleTreatments(x, rules)
  ret <- list(x = x, y = y, A = A, oracle = oracle, treatment.effects = treatment.effects)
  class(ret) <- "sim.owl.data"
  ret
}

simulateOwlData2 <- function(n, p, rules, true.beta, interaction, sd.x = 1, sd.y = 1,
                             outcome.type = c("numeric", "binary", "misspecified.binary"), 
                             num.factors = 0L, factor.levels = rep(3, num.factors), 
                             misspecified = FALSE) {
  
  outcome.type <- match.arg(outcome.type)
  if (num.factors > 0) stopifnot(all(factor.levels > 1))
  num.contin <- p - num.factors
  num.factors <- as.integer(num.factors)
  n <- as.integer(n)
  p <- as.integer(p)
  x <- matrix(rnorm(n * num.contin, sd = sd.x), n, num.contin)
  colnames(x) <- paste("C.X", 1:num.contin, sep = "")
  
  if (num.factors > 0){
    num.binary <- sum(factor.levels == 2)
    
    #x <- matrix(rbinom(n * p, 1, 0.05), n, p)
    x.factors <- sapply(factor.levels, function(x) as.factor(sample.int(x, n, replace = TRUE)))
    colnames(x.factors) <- paste("F.X", 1:num.factors, sep = "")
    colnames(x.factors)[which(factor.levels == 2)] <- paste(paste("Bin.X", 1:num.binary, sep = ""), ".", sep = "")
    colnames(x.factors)[which(factor.levels != 2)] <- paste(paste("F.X", 1:(num.factors - num.binary), sep = ""), ".", sep = "")
    
    x <- data.frame(x, x.factors)
  } else {
    x <- data.frame(x)
  }
  
  
  if (interaction) {
    x <- model.matrix(~ . + .*., data = x)[,-1]
    #x <- modelMatrixNeg(x, interaction = TRUE)
    #non.int <- colnames(x)[-grep(":", colnames(x))]
    all.uniques <- gsub("[\\.]([0-9]+)", "", colnames(x))
    contins <- grep("^(C|Bin)\\.X[0-9]+[^F]*$[^:]*$", all.uniques, perl = TRUE)
    #uniques <- gsub("[\\.]([0-9]+)", "", non.int)
    #factors.non.int <- colnames(x)[-grep(":|[C]|[Bin]", colnames(x))]
    grouping <- rep(NA, ncol(x))
    uni <- unique(all.uniques[-contins])
    grouping[-contins] <- match(all.uniques[-contins], uni)
  } else {
    x <- model.matrix(~ ., data = x)[,-1]
    #x <- modelMatrixNeg(x, interaction = FALSE)
    all.uniques <- gsub("[\\.]([0-9]+)", "", colnames(x))
    contins <- grep("^(C|Bin)\\.X[0-9]+", all.uniques, perl = TRUE)
    grouping <- rep(NA, ncol(x))
    uni <- unique(all.uniques[-contins])
    grouping[-contins] <- match(all.uniques[-contins], uni)
  }
  
  A <- as.factor(sample(3, n, replace = TRUE))
  
  treatment.effects <- genTreatmentEffects2(x, A, rules)
  
  if (misspecified) {
    y.main.eff <- rnorm(n * length(true.beta), ncol = length(true.beta)) %*% true.beta
  } else {
    y.main.eff <- x[,3:(3 + length(true.beta) - 1)] %*% true.beta
  }
  
  if (outcome.type == "binary") {
    log.p.ratio <- rowSums(treatment.effects) + y.main.eff
    prob.y.1 <- 1 / (1 + exp(-log.p.ratio))
    y <- rbinom(n, 1, prob = prob.y.1)
  } else if (outcome.type == "numeric") {
    y <- rnorm(n, sd = sd.y) + rowSums(treatment.effects) + y.main.eff
    y <- y + abs(min(y))
  } else if (outcome.type == "misspecified.binary"){
    exp.p.ratio <- rowSums(treatment.effects) + y.main.eff + 1
    prob.y.1 <- log(exp.p.ratio + sign(min(exp.p.ratio)) * min(exp.p.ratio) + 1) / (1 + log(exp.p.ratio + sign(min(exp.p.ratio)) * min(exp.p.ratio) + 1))
    y <- rbinom(n, 1, prob = prob.y.1)
  }
  
  oracle <- genOracleTreatments2(x, rules)
  ret <- list(x = x, y = y, A = A, oracle = oracle, treatment.effects = treatment.effects, grouping = grouping)
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