
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
