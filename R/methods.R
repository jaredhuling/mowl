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
  
  n.trt <- nrow(x$d.vals)
  n.lam <- ncol(x$d.vals)
  d.vals <- numeric(length = (n.trt) * n.lam)
  for (i in 1:(nrow(x$d.vals))) {
    d.vals[((i-1) * n.lam + 1):(i * n.lam)] <- x$d.vals[i,]
  }
  
  dfs <- if(inherits(x$model, "msgl")) {df.msgl(x)} else {x$model$df}
  
  dvaldat <- data.frame(dvals = d.vals, 
                        treatment = rep(as.character(1:4), each = ncol(x$d.vals)),
                        df = rep(dfs, nrow(x$d.vals)),
                        lambda = rep(x$model$lambda, nrow(x$d.vals)))
  
  dfdat <- data.frame(df = rep(dfs, nrow(x$d.vals)),
                      lambda = rep(x$model$lambda, nrow(x$d.vals)))
  print(x$class.lambda)
  class.lam <- -1 * x$class.lambda
  value.lam <- -1 * x$value.lambda
  print(value.lam)
  aic.lam <- -1 * x$aic.lambda
  class.text.y.p2 <- (x$model$df[x$class.lambda.idx] + max(x$model$df)) / 2.5
  value.text.y.p2 <- (x$model$df[x$value.lambda.idx] + max(x$model$df)) / 2.5
  aic.text.y.p2 <- (x$model$df[x$aic.lambda.idx] + max(x$model$df)) / 2.5
  value.text <- if(class.lam == value.lam) {""} else ("\n Value Func. Selection")
  
  p1 <- ggplot(aes(x = -lambda, y = dvals, color = treatment), data = dvaldat) + geom_line(size=1.2) + 
    theme_bw() + theme(legend.position = "bottom") + geom_vline(xintercept = -class.lam) +
    geom_vline(xintercept = -aic.lam)
  

  
  p2 <- ggplot(aes(x = -lambda, y = df), data = dfdat) + geom_line(size=1.2) + theme_bw() +
    geom_vline(xintercept = class.lam) +
    geom_text(aes(x = class.lam, label="\n Misclassification Selection", y = class.text.y.p2), 
              colour="blue", angle=90, text=element_text(size=11)) +
    geom_vline(xintercept = aic.lam) +
    geom_text(aes(x = aic.lam, label="\n AIC Selection", y = aic.text.y.p2), 
              colour="blue", angle=90, text=element_text(size=11)) + 
    geom_vline(xintercept = value.lam) +
    geom_text(aes(x = value.lam, label=value.text, y = value.text.y.p2), 
              colour="blue", angle=90, text=element_text(size=11))
  print(grid.draw(rbind(ggplotGrob(p1), ggplotGrob(p2), size = "last")))
}

