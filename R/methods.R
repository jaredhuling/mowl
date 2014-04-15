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
  
  if (!is.null(obj$max.pct.correct)) {
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
  
  dvals <- data.frame(array(0, dim = c(4, (length(obj$d.optimal)+1))))
  rownames(dvals) <- c("Optimal", "Class", "Value", "AIC")
  colnames(dvals) <- c("Lambda", paste("D", 1:length(obj$d.optimal), sep = ""))
  dvals[1, 2:(length(obj$d.optimal)+1)] <- obj$d.optimal
  dvals[2, 2:(length(obj$d.optimal)+1)] <- obj$d.class
  dvals[3, 2:(length(obj$d.optimal)+1)] <- obj$d.value
  dvals[4, 2:(length(obj$d.optimal)+1)] <- obj$d.aic
  dvals[1, 1] <- obj$optimal.d.lambda
  dvals[2,1] <- obj$class.lambda
  dvals[3,1] <- obj$value.lambda
  dvals[4,1] <- obj$aic.lambda
  print(dvals)
  
}


plot.owlfit <- function(x) {
  #K <- nrow(obj$d.vals)
  #nonan <- obj$d.vals; nonan[is.nan(nonan)] <- mean(nonan[!is.nan(nonan)])
  #ylims <- c(min(na.omit(obj$d.vals))-0.25 * sd(nonan[!is.nan(nonan)]), max(na.omit(obj$d.vals)+0.25 * sd(nonan[!is.nan(nonan)])))
  #cols <- c("red", "green", "yellow", "purple", "orange", colors()[sample.int(400, K)])
  #plot(x = 1:length(obj$model$lambda), y = obj$d.vals[1,], type = "l", col = "blue", ylim = ylims)
  #for (i in 2:K) {
  #  lines(x = 1:length(obj$model$lambda), y = obj$d.vals[i,], col = cols[i-1])
  #}
  #print(ylims)
  require(grid)
  require(gridExtra)
  
  n.trt <- nrow(obj$d.vals)
  n.lam <- ncol(obj$d.vals)
  d.vals <- numeric(length = (n.trt) * n.lam)
  for (i in 1:(nrow(obj$d.vals))) {
    d.vals[((i-1) * n.lam + 1):(i * n.lam)] <- obj$d.vals[i,]
  }
  
  dfs <- if(inherits(obj$model, "msgl")) {df.msgl(obj)} else {obj$model$df}
  
  dvaldat <- data.frame(dvals = d.vals, 
                        treatment = rep(as.character(1:4), each = ncol(obj$d.vals)),
                        df = rep(dfs, nrow(obj$d.vals)),
                        lambda = rep(obj$model$lambda, nrow(obj$d.vals)))
  
  dfdat <- data.frame(df = rep(dfs, nrow(obj$d.vals)),
                      lambda = rep(obj$model$lambda, nrow(obj$d.vals)))
  
  p1 <- ggplot(aes(x = -lambda, y = dvals, color = treatment), data = dvaldat) + geom_line(size=1.2) + 
    theme_bw() + theme(legend.position = "bottom") + geom_vline(xintercept = -obj$class.lambda) +
    geom_vline(xintercept = -obj$aic.lambda)
  
  p2 <- ggplot(aes(x = -lambda, y = df), data = dfdat) + geom_line(size=1.2) + theme_bw() +
    geom_vline(xintercept = -obj$class.lambda) +
    geom_text(aes(x = -obj$class.lambda, label="\n Misclassification Selection", y = (obj$model$df[obj$class.lambda.idx] + max(obj$model$df)) / 2.5), 
              colour="blue", angle=90, text=element_text(size=11)) +
    geom_vline(xintercept = -obj$aic.lambda) +
    geom_text(aes(x = -obj$aic.lambda, label="\n AIC Selection", y = (obj$model$df[obj$aic.lambda.idx] + max(obj$model$df)) / 2.5), 
              colour="blue", angle=90, text=element_text(size=11))
  print(grid.draw(rbind(ggplotGrob(p1), ggplotGrob(p2), size = "last")))
  list(p1 = p1, p2 = p2)
}

