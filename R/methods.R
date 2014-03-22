predict.owlfit <- function(object, type.measure = c("class", "value", "aic", "deviance", "optimal"), ...) {
  type.measure <- match.arg(type.measure)
  lam <- switch(type.measure,
                class = object$class.lambda,
                value = object$value.lambda,
                aic = object$aic.lambda,
                optimal = object$optimal.lambda)
  predict(object$model, s = lam, ...)
}

print.owlfit <- function(obj) {
  cat("\nCall: ", deparse(obj$call), "\n\n")
  
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
  cat("Optimal d \n\n")
  print(obj$d.optimal)
  cat("Misclassification Criterion d \n\n")
  print(obj$d.class)
  cat("Value Function Criterion d \n\n")
  print(obj$d.value)
  cat("AIC Criterion d \n\n")
  print(obj$d.aic)
}