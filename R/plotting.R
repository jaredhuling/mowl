
genContourData <- function(x, rules) {
  contour.data <- expand.grid(x1 = seq(min(x[,1]), max(x[,1]), length = 300),
                              x2 = seq(min(x[,2]), max(x[,2]), length = 300))
  contour.y <- genOracleTreatments(contour.data, rules)
  data.frame(contour.data, z = as.numeric(contour.y))
}


plotTreatmentRules <- function(x, predictions, rules) {
  require(ggplot2)
  contour.data <- genContourData(x, rules)
  dat2plot <- data.frame(x1 = x[,1], x2 = x[,2], y = predictions); names(dat2plot) <- c("x1", "x2", "y")
  
  pl <- ggplot(aes( x1 , x2 ), data = dat2plot ) + 
    geom_point( aes(colour = factor(y), shape = factor(y), alpha = 0.75), size = 10) + 
    theme_bw()+ scale_shape_manual(name = "Treatment", values = as.character(1:3)) +
    scale_colour_discrete(name = "Treatment") + theme(legend.text = element_blank())
  pl + stat_contour(aes( x1, x2, z = z), data = contour.data, size = 1) + scale_alpha(guide = 'none')
  
}

plotPatientEffects <- function(obj, x, patient.ind, lam.ind, return.all = FALSE) {
  require(ggplot2)
  require(grid)
  require(gridExtra)
  
  effects <- trts <- NULL
  if (inherits(obj$model, "glmnet")) {
    for (i in 1:length(obj$model$beta)) {
      c.betas <- obj$model$beta[[i]][,lam.ind]
      nz.ind <- which(c.betas != 0)
      c.betas.nz <- c.betas[nz.ind]
      col.ind <- match(names(c.betas.nz), colnames(x))
      effects.tmp <- c(obj$model$a0[i, lam.ind], x[patient.ind, col.ind] * c.betas.nz)
      effects.tmp <- effects.tmp[order(abs(effects.tmp), decreasing = TRUE)]
      effects <- c(effects, effects.tmp)
      trts <- c(trts, rep(names(obj$model$beta)[i], length(effects.tmp)))
    }
  } else {
    full.beta <- as(obj$model$beta[[lam.ind]], "matrix")
    trt.names <- dimnames(full.beta)[[1]]
    for (i in 1:nrow(obj$model$beta[[1]])) {
      c.betas <- full.beta[i,-1]
      nz.ind <- which(c.betas != 0)
      c.betas.nz <- c.betas[nz.ind]
      col.ind <- match(names(c.betas.nz), colnames(x))
      effects.tmp <- c(full.beta[i,1], x[patient.ind, col.ind] * c.betas.nz)
      effects.tmp <- effects.tmp[order(abs(effects.tmp), decreasing = TRUE)]
      effects <- c(effects, effects.tmp)
      trts <- c(trts, rep(trt.names[i], length(effects.tmp)))
    }
  }
  
  dat2plot <- data.frame(Treatment = trts, effects = effects)
  dat2plotpos <- subset(dat2plot,effects >= 0)
  dat2plotneg <- subset(dat2plot,effects < 0)
  total.effects <- tapply(dat2plot$effects, as.factor(dat2plot$Treatment), FUN = sum)
  pl <- ggplot() + 
    geom_bar(data = dat2plotpos, aes(x=Treatment, y=effects, 
                                     fill=factor(effects)),stat = "identity") +
    geom_bar(data = dat2plotneg, aes(x=Treatment, y=effects, 
                                     fill=factor(effects)),stat = "identity") +
    theme_bw() + theme(legend.position="none",axis.title.x=element_blank(),
                       plot.margin = grid::unit(c(1, 1, 0, 0.5), "cm"))
  dat2plot2 <- data.frame(Treatment = names(total.effects), effects = 
                            total.effects)
  pl2 <- ggplot() + 
    geom_bar(data = dat2plot2, aes(x=Treatment, y=effects),stat = "identity") +
    theme_bw() + theme(legend.position="none", plot.margin = grid::unit(c(0, 0.5, 0.5, 0.15), "cm"))
  #print(pl)
  gp1<- ggplot_gtable(ggplot_build(pl))
  gp2<- ggplot_gtable(ggplot_build(pl2))
  print(grid.draw(rbind(ggplotGrob(pl), ggplotGrob(pl2), size="last")))
  if (return.all) {
    list(plot = pl, dat2plot = dat2plot)
  } else {
    list(plot = pl, dat2plot = dat2plot)
  }
}


plotTreatmentEffects <- function(obj, lam.ind = NULL) {
  if (is.null(lam.ind)){lam.ind <- obj$value.lambda.idx}
  n.trts <- nrow(obj$d.vals)
  dat.list <- vector(mode = "list", length = n.trts)
  plot.obj <- ggplot(mapping = aes(x=variable, y=value)) + theme_bw() +
    theme(axis.text.x = element_text(angle=90, hjust=1, size = 14))
  
  if (inherits(obj$model, "msgl")) {
    full.beta <- as(obj$model$beta[[lam.ind]], "matrix")[,-1]
  }
  
  for (i in 1:n.trts) {
    if (inherits(obj$model, "msgl")) {
      c.betas <- full.beta[i, ]
    } else {
      c.betas <- obj$model$beta[[i]][,lam.ind]
    }
    c.betas.nz <- c.betas[c.betas !=0]
    t1 <- c.betas.nz[order(abs(c.betas.nz), decreasing = TRUE)]
    
    ## Plot Treatment 1 effects
    dat2plott1 <- data.frame(variable = names(t1), value = t1, treatment = rep(as.character(i), length(t1)))
    
    dat2plott1 <- transform(dat2plott1, variable=reorder(variable, -abs(value)) ) 
    
    dat2plott1$variable <- factor(dat2plott1$variable, levels=unique(as.character(dat2plott1$variable)), ordered = TRUE )
    dat.list[[i]] <- dat2plott1
    plot.obj <- plot.obj + geom_point(data=dat2plott1, size = 4)
  }
  plot.obj <- plot.obj + facet_wrap(~ treatment, scale = "free_x", ncol = 2) 
  
  print(plot.obj)
}



