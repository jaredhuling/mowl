predict.owlfit <- function(object, type.measure = c("class", "value", "aic", "deviance", "optimal"), ...) {
  type.measure <- match.arg(type.measure)
  lam <- switch(type.measure,
                class = object$class.lambda,
                value = object$value.lambda,
                aic = object$aic.lambda,
                optimal = object$optimal.lambda)
  if (inherits(object$model, "msgl")) {
    predict(object$model, ...)$classes[,which(object$model$lambda == lam)]
  } else {
    predict(object$model, s = lam, ...)
  }
}

print.owlfit <- function(obj) {
  cat("\nCall: ", deparse(obj$call), "\n\n")
  
  if (!is.null(max.pct.correct)) {
    lams <- data.frame(array(0, dim = c(4, 2)))
    rownames(lams) <- c("Optimal", "Class", "Value", "AIC")
    colnames(lams) <- c("Lambda", "Pct Correct")
    lams[1,1] <- obj$optimal.lambda
    lams[1,2] <- obj$max.pct.correct
    lams[2,1] <- obj$class.lambda
    lams[2,2] <- obj$class.pct.correct
    lams[3,1] <- obj$value.lambda
    lams[3,2] <- obj$value.pct.correct
    lams[4,1] <- obj$aic.lambda
    lams[4,2] <- obj$aic.pct.correct
    print(lams)
  }
  
  cat("\n")
  
  dvals <- data.frame(array(0, dim = c(4, 4)))
  rownames(dvals) <- c("Optimal", "Class", "Value", "AIC")
  colnames(dvals) <- c("Lambda", "D1", "D2", "D3")
  dvals[1, 2:4] <- obj$d.optimal
  dvals[2, 2:4] <- obj$d.class
  dvals[3, 2:4] <- obj$d.value
  dvals[4, 2:4] <- obj$d.aic
  dvals[1, 1] <- obj$optimal.d.lambda
  dvals[2,1] <- obj$class.lambda
  dvals[3,1] <- obj$value.lambda
  dvals[4,1] <- obj$aic.lambda
  print(dvals)
  
}


plot.owlfit <- function(obj) {
  K <- nrow(obj$d.vals)
  nonan <- obj$d.vals; nonan[is.nan(nonan)] <- mean(nonan[!is.nan(nonan)])
  ylims <- c(min(na.omit(obj$d.vals))-0.25 * sd(nonan[!is.nan(nonan)]), max(na.omit(obj$d.vals)+0.25 * sd(nonan[!is.nan(nonan)])))
  cols <- c("red", "green", "yellow", "purple", "orange", colors()[sample.int(400, K)])
  plot(x = 1:length(obj$model$lambda), y = obj$d.vals[1,], type = "l", col = "blue", ylim = ylims)
  for (i in 2:K) {
    lines(x = 1:length(obj$model$lambda), y = obj$d.vals[i,], col = cols[i-1])
  }
  print(ylims)
}

