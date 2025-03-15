xMoments <- function(x, p) {
  n <- nrow(x)
  m <- ncol(x)
  mu1 <- array(0, m)
  mu2 <- array(0, c(m, m))
  mu3 <- array(0, c(m, m, m))
  mu4 <- array(0, c(m, m, m, m))
  for (i in 1:n) {
    mu1 <- mu1 + p[i] * x[i, ]
    mu2 <- mu2 + p[i] * outer(x[i, ], x[i, ])
    mu3 <- mu3 + p[i] * outer(outer(x[i, ], x[i, ]), x[i, ])
    mu4 <- mu4 + p[i] * outer(outer(outer(x[i, ], x[i, ]), x[i, ]), x[i, ])
  }
  return(list(
    mu1 = mu1,
    mu2 = mu2,
    mu3 = mu3,
    mu4 = mu4
  ))
}

multinomialMoments <- function(n, p) {
  m <- length(p)
  mu1 <- array(0, m)
  mu2 <- array(0, c(m, m))
  mu3 <- array(0, c(m, m, m))
  mu4 <- array(0, c(m, m, m, m))
  ff <- fallingFactorials(n, 4)
  for (i in 1:m) {
    mu1[i] <- n * p[i]
    for (j in 1:m) {
      mu2[i, j] <- ff[2] * p[i] * p[j] + n * kroneckerDelta(c(i, j)) * p[i]
      for (k in 1:m) {
        mu3[i, j, k] <- ff[3] * p[i] * p[j] * p[k]
        mu3[i, j, k] <- mu3[i, j, k] + ff[2] * kroneckerDelta(c(i, k)) * p[j] * p[k]
        mu3[i, j, k] <- mu3[i, j, k] + ff[2] * kroneckerDelta(c(j, k)) * p[i] * p[k]
        mu3[i, j, k] <- mu3[i, j, k] + ff[2] * kroneckerDelta(c(i, j)) * p[j] * p[k]
        mu3[i, j, k] <- mu3[i, j, k] + n * kroneckerDelta(c(i, j, k)) * p[i]
        for (l in 1:m) {
          mu4[i, j, k, l] <- ff[4] * p[i] * p[j] * p[k] * p[l]
          mu4[i, j, k, l] <- mu4[i, j, k, l] + ff[3] * kroneckerDelta(c(k, l)) * p[i] * p[j] * p[l]
          mu4[i, j, k, l] <- mu4[i, j, k, l] + ff[3] * kroneckerDelta(c(j, l)) * p[i] * p[k] * p[l]
          mu4[i, j, k, l] <- mu4[i, j, k, l] + ff[3] * kroneckerDelta(c(i, l)) * p[j] * p[k] * p[l]
          mu4[i, j, k, l] <- mu4[i, j, k, l] + ff[3] * kroneckerDelta(c(i, k)) * p[j] * p[k] * p[l]
          mu4[i, j, k, l] <- mu4[i, j, k, l] + ff[3] * kroneckerDelta(c(j, k)) * p[i] * p[k] * p[l]
          mu4[i, j, k, l] <- mu4[i, j, k, l] + ff[3] * kroneckerDelta(c(i, j)) * p[i] * p[k] * p[l]
          mu4[i, j, k, l] <- mu4[i, j, k, l] + ff[2] * kroneckerDelta(c(i, k, l)) * p[j] * p[l]
          mu4[i, j, k, l] <- mu4[i, j, k, l] + ff[2] * kroneckerDelta(c(j, k, l)) * p[i] * p[l]
          mu4[i, j, k, l] <- mu4[i, j, k, l] + ff[2] * kroneckerDelta(c(i, j, l)) * p[k] * p[l]
          mu4[i, j, k, l] <- mu4[i, j, k, l] + ff[2] * kroneckerDelta(c(i, j, k)) * p[k] * p[l]
          mu4[i, j, k, l] <- mu4[i, j, k, l] + ff[2] * kroneckerDelta(c(i, k)) * kroneckerDelta(c(j, l)) * p[k] * p[l]
          mu4[i, j, k, l] <- mu4[i, j, k, l] + ff[2] * kroneckerDelta(c(j, k)) * kroneckerDelta(c(i, l)) * p[k] * p[l]
          mu4[i, j, k, l] <- mu4[i, j, k, l] + ff[2] * kroneckerDelta(c(i, j)) * kroneckerDelta(c(k, l)) * p[j] * p[l]
          mu4[i, j, k, l] <- mu4[i, j, k, l] + n * kroneckerDelta(c(i, j, k, l)) * p[l]
        }
      }
    }
  }
  return(list(
    mu1 = mu1,
    mu2 = mu2,
    mu3 = mu3,
    mu4 = mu4
  ))
}

kroneckerDelta <- function(i) {
  return(ifelse(all(i == i[1]), 1, 0))
}

fallingFactorials <- function(n, k) {
  return(factorial(n) / factorial(n - 1:k))
}
