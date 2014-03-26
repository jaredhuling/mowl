
value.func <- function(A, predicted, y) {
  A.char <- levels(A)[A]
  ind <- 1 * (A.char == predicted)
  sum(3 * y * ind) / length(A)
}


logLikNull <- function(y, A, weights) {
  require(mlogit)
  mdat <- data.frame(A = A, weights = weights)
  rownames(mdat) <- 1:length(weights)
  mdat2 <- mlogit.data(mdat, choice = "A", shape = "wide")
  mfit <- mlogit(A ~ 1, weights = weights, data = mdat2)
  logLik(mfit)[1]
}

logLikGlmnet <- function(y, obj, A, weights) {
  logLikSat <- (obj$nulldev / 2) + logLikNull(y, A, weights)
  logLikSat - deviance(obj) / 2
}

aic <- function(obj, y, A, weights) {
  stopifnot(inherits(obj, "glmnet"))
  -2 * logLikGlmnet(y, obj, A, weights) + 2 * obj$df
}

df.msgl <- function(obj) {
  betas <- array(0, dim = length(obj$beta))
  for (i in 1:length(obj$beta)) {
    betas[i] <- sum(colSums(obj$beta[[i]]) != 0)
  }
  betas
}

aic.msgl <- function(obj) {
  2 * (df.msgl(obj) + obj$loss * length(obj$classes.true))
}

bic <- function(obj, y, A, weights) {
  stopifnot(inherits(obj, "glmnet"))
  -2 * logLikGlmnet(y, obj, A, weights) + obj$df * (log(obj$nobs) + log(2 * pi))
}


predTreatment <- function(obj, x) {
  if (inherits(obj, "cv.glmnet")) {
    obj <- obj$glmnet.fit
  }
  coef.vals <- array(0, dim = c(length(obj$beta), nrow(x), ncol(obj$beta[[1]])))
  for (i in 1:length(obj$beta)) {
    coef.vals[i, , ] <- drop(cbind(rep(1, nrow(x)), x) %*% rbind(obj$a0[i,], as(obj$beta[[i]], "matrix")))
  }
  apply(coef.vals, 2:3, function(x) as.character(which.max(x)))
}
