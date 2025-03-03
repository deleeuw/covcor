library(combinat)

make4 <- function(m, p = rep(1, 4) / 4) {
  matje <- NULL
  pitje <- NULL
  for (i in 0:m) {
    for (j in 0:m) {
      for (k in 0:m) {
        for (l in 0:m) {
          vec <- c(i, j, k, l)
          if (sum(vec) == m) {
            matje <- rbind(matje, vec, deparse.level = 0)
            pitje <- c(pitje, dmnom(vec, prob = p))
          }
        }
      }
    }
  }
  return(list(matje = matje, pitje = pitje))
}

triple <- function(m, p, i, j, k, l) {
  z <- make4(m, p)
  xx <- z$matje
  pp <- z$pitje
  mu <- colSums(pp * xx)
  for (nu in 1:nrow(xx)) {
    xx[nu, ] <- xx[nu, ] - mu
  }
  np <- length(pp)
  sijkl <- sij <- skl <- 0.0
  for (nu in 1:np) {
    sijkl <- sijkl + pp[i] * xx[nu, i] * xx[nu, j] * xx[nu, k] * xx[nu, l]
    sij <- sij + pp[nu] * xx[nu, i] * xx[nu, j]
    skl <- skl + pp[nu] * xx[nu, k] * xx[nu, l]
  }
  return(c(sijkl, sij, skl, sijkl - sij * skl))
}