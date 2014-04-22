
groups <- function(x) {
  attr(x, "groups")
}

thresholdModel <- function(obj, threshold = 0) {
  owl.obj <- obj
  if (inherits(obj, "owlfit")) {
    owlfit <- TRUE
  } else if (inherits(obj, "msgl") | inherits(obj, "glmnet")) {
    owlfit <- FALSE
  } else {
    stop("need either owlfit, msgl, or glmnet object")
  }
  
  if (owlfit) {
    nbeta <- length(owl.obj$model$beta)
    for (i in 1:nbeta) {
      owl.obj$model$beta[[i]][abs(owl.obj$model$beta[[i]]) < threshold] <- 0
    }
  } else {
    nbeta <- length(owl.obj$beta)
    for (i in 1:nbeta) {
      owl.obj$beta[[i]][abs(owl.obj$beta[[i]]) < threshold] <- 0
    }
  }

  owl.obj
}


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

print.owlfit <- function(obj, xtable = FALSE,  ...) {
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
  if (xtable) {
    require(xtable)
    xt.res <- xtable(dvals, digits = 3)
    align(xt.res) <- paste("r", paste(rep("l", (length(obj$d.optimal) + 1)), collapse = ""), sep = "|")
    print(xt.res, floating = FALSE)
  }
  
}


plot.owlfit <- function(x, log.scale = FALSE) {
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
  d.vals <- numeric(length = (n.trt * n.lam))
  for (i in 1:(nrow(x$d.vals))) {
    d.vals[((i-1) * n.lam + 1):(i * n.lam)] <- x$d.vals[i,]
  }
  
  dfs <- if(inherits(x$model, "msgl")) {df.msgl(x$model)} else {x$model$df}
  
  if (log.scale) {
    class.lam <- -1 * log(x$class.lambda)
    value.lam <- -1 * log(x$value.lambda)
    aic.lam <- -1 * log(x$aic.lambda)
    lams <- rep(log(x$model$lambda), n.trt)
  } else {
    class.lam <- -1 * x$class.lambda
    value.lam <- -1 * x$value.lambda
    aic.lam <- -1 * x$aic.lambda
    lams <- rep(x$model$lambda, n.trt)
  }
  
  dvaldat <- data.frame(dvals = d.vals, 
                        treatment = rep(as.character(1:4), each = n.lam),
                        df = rep(dfs, n.trt),
                        lambda = lams)
  
  dfdat <- data.frame(df = rep(dfs, n.trt),
                      lambda = lams)

  class.text.y.p2 <- (dfs[x$class.lambda.idx] + max(dfs)) / 2.5
  value.text.y.p2 <- (dfs[x$value.lambda.idx] + max(dfs)) / 2.5
  aic.text.y.p2 <- (dfs[x$aic.lambda.idx] + max(dfs)) / 2.5
  value.text <- if(class.lam == value.lam) {""} else ("\n Value Func. Selection")
  
  vline.dat <- data.frame(value.lam = value.lam, class.lam = class.lam,
                          aic.lam = aic.lam, value.text = value.text,
                          class.text.y.p2 = class.text.y.p2,
                          value.text.y.p2 = value.text.y.p2,
                          aic.text.y.p2 = aic.text.y.p2)
  
  x.label <- if (log.scale) {"-log(lambda)"} else {"-lambda"}
  
  p1 <- ggplot(aes(x = -lambda, y = dvals, color = treatment), data = dvaldat) + geom_line(size=1.2) + 
    theme_bw() + theme(legend.position = "bottom") + xlab(x.label) +
    geom_vline(aes(xintercept = class.lam), data = vline.dat) +
    geom_vline(aes(xintercept = value.lam), data = vline.dat) +
    geom_vline(aes(xintercept = aic.lam), data = vline.dat)
  

  
  p2 <- ggplot(aes(x = -lambda, y = df), data = dfdat) + geom_line(size=1.2) + theme_bw() +
    geom_vline(aes(xintercept = class.lam), data = vline.dat) + xlab(x.label) +
    geom_text(aes(x = class.lam, label="\n Misclassification Selection", y = class.text.y.p2),
              data = vline.dat,
              colour="blue", angle=90, text=element_text(size=11)) +
    geom_vline(aes(xintercept = aic.lam), data = vline.dat) +
    geom_text(aes(x = aic.lam, label="\n AIC Selection", y = aic.text.y.p2), data = vline.dat,
              colour="blue", angle=90, text=element_text(size=11)) + 
    geom_vline(aes(xintercept = value.lam), data = vline.dat) +
    geom_text(aes(x = value.lam, label=value.text, y = value.text.y.p2), data = vline.dat,
              colour="blue", angle=90, text=element_text(size=11))
  print(grid.draw(rbind(ggplotGrob(p1), ggplotGrob(p2), size = "last")))
}

