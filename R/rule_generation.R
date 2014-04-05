returnRule <- function(coefs) {
  #returns function that evaluates treatment rule
  if (length(coefs) == 4) {
    function(x1 = NULL, x2 = NULL, solve = FALSE) {
      c0 <- coefs[1]
      c1 <- coefs[2]
      c2 <- coefs[3]
      c3 <- coefs[4]
      #print(coefs)
      if (solve) {
        c(c0, c1, c2)
      } else {
        c0 + c1 * x1 + c2 * x2 + c3 * x1 * x2
      }
    }
  } else if (length(coefs) == 3) {
    function(x1 = NULL, x2 = NULL, solve = FALSE) {
      c0 <- coefs[1]
      c1 <- coefs[2]
      c2 <- coefs[3]
      #print(coefs)
      if (solve) {
        c(c0, c1, c2)
      } else {
        c0 + c1 * x1 + c2 * x2
      }
    }
  } else {stop("Wrong number of rule coefficients")}
}

returnRule2 <- function(lst) {
  #returns function that evaluates treatment rule
  print(lst)
  function(x) {
    
    print(lst)
    ret <- x[,lst$var.idx] %*% lst$coefs
    int <- rowSums(apply(lst$int.idx, 2, function(idx) lst$int.coefs * x[,idx[1]] * x[,idx[2]]))
    ret + int
    
  }

}


# input rules must be a list of length equal to the number of treatments.
# each element must be a list with: 
#     'coefs' for coefficient values of treatment effect
#     'var.idx' for index of variables corresponding to 'coefs'
#     'int.coefs' for coefficient values for interaction effects on treatment
#     'int.idx' is a matrix with number of columns equal to length of 'int.coefs'
#               and for each column, the first value is the index corresponding to
#               the variable index for first part of interaction and second value for 
#               the second part of interaction
genTreatmentRules2 <- function(lst) {
  stopifnot(class(lst) == "list")
  rules <- vector(mode = "list", length = length(lst))
  for (i in 1:length(lst)) {
    rules[[i]] <- returnRule2(lst[[i]])
  }
  names(rules) <- paste("rule", 0:(length(lst)-1), sep = "")
  class(rules) <- unique(c(class(rules), "treatment.rules"))
  rules
}

genTreatmentRules <- function(lst) {
  stopifnot(class(lst) == "list")
  rules <- vector(mode = "list", length = length(lst))
  for (i in 1:length(lst)) {
    rules[[i]] <- returnRule(lst[[i]])
  }
  names(rules) <- paste("rule", 0:(length(lst)-1), sep = "")
  class(rules) <- unique(c(class(rules), "treatment.rules"))
  rules
}

solve.treatment.rules <- function(rules) {
  A <- array(0, dim = c(length(rules), 3))
  B <- array(0, dim = c(length(rules), 3))
  k <- 0
  for (i in 1:(length(rules) - 1)) {
    for (j in (i + 1):length(rules)) {
      k <- k + 1
      A[k, ] <- rules[[i]](solve = TRUE)
      B[k, ] <- rules[[j]](solve = TRUE)
    }
  }
  list(A = A, B = B)
}


genOracleTreatments <- function(x, rules) {
  mat <- array(0, dim = c(nrow(x), length(rules)))
  for (i in 1:ncol(mat)) {
    mat[,i] <- rules[[i]](x[,1], x[,2])
  }
  apply(mat, 1, function(xx) as.character(which.max(xx)))
}

genOracleTreatments2 <- function(x, rules) {
  mat <- array(0, dim = c(nrow(x), length(rules)))
  for (i in 1:ncol(mat)) {
    mat[,i] <- rules[[i]](x)
  }
  apply(mat, 1, function(xx) as.character(which.max(xx)))
}

genTreatmentEffects <- function(x, A, rules) {
  stopifnot(inherits(A, "factor"))
  effects <- array(0, dim = c(length(A), length(rules)))
  A.mm <- model.matrix( ~ A - 1)
  for (i in 1:length(rules)) {
    effects[,i] <- rules[[i]](x[,1], x[,2]) * A.mm[,i]
  }
  effects
}

genTreatmentEffects2 <- function(x, A, rules) {
  stopifnot(inherits(A, "factor"))
  effects <- array(0, dim = c(length(A), length(rules)))
  A.mm <- model.matrix( ~ A - 1)
  for (i in 1:length(rules)) {
    effects[,i] <- rules[[i]](x) * A.mm[,i]
  }
  effects
}

