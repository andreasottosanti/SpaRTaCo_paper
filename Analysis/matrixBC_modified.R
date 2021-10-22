library(sparseBC)

matrixBC2 <- function (x, k, r, lambda, alpha, beta, nstart = 20, Cs.init = NULL,
                       Ds.init = NULL, max.iter = 50, threshold = 1e-04, Sigma.init = NULL,
                       Delta.init = NULL, center = TRUE){

    MatrixObjective2 <- function (x, mus, Cs, Ds, Sigma, Delta, lambda = 0, alpha = 0,
                                  beta = 0)
    {
      tempSigma <- 0
      tempSigma2 <- 0
      tempDelta <- 0
      tempDelta2 <- 0
      tempsum <- 0
      tempmus <- 0
      for (i in 1:max(Cs)) {
        temp1 <- determinant(solve(Sigma[Cs == i, Cs == i]), logarithm = T)$modulus
        tempSigma <- temp1 + tempSigma
        temp11 <- solve(Sigma[Cs == i, Cs == i])
        tempSigma2 <- sum(abs(temp11)) + tempSigma2
      }
      for (j in 1:max(Ds)) {
        temp2 <- determinant(solve(Delta[Ds == j, Ds == j]), logarithm = T)$modulus
        tempDelta <- temp2 + tempDelta
        temp21 <- solve(Delta[Ds == j, Ds == j])
        tempDelta2 <- sum(abs(temp21)) + tempDelta2
      }
      for (i in 1:max(Cs)) {
        for (j in 1:max(Ds)) {
          temp3 <- solve(Sigma[Cs == i, Cs == i]) %*% (matrix(x[Cs ==
                                                                  i, Ds == j], nrow = sum(Cs == i), ncol = sum(Ds ==
                                                                                                                 j)) - mus[Cs == i, Ds == j]) %*% solve(Delta[Ds ==
                                                                                                                                                                j, Ds == j]) %*% t(matrix(x[Cs == i, Ds == j],
                                                                                                                                                                                          nrow = sum(Cs == i), ncol = sum(Ds == j)) - mus[Cs ==
                                                                                                                                                                                                                                            i, Ds == j])
          tempsum <- sum(diag(temp3)) + tempsum
          temp4 <- abs(mus[Cs == i, Ds == j][1])
          tempmus <- temp4 + tempmus
        }
      }
      obj <- -(ncol(x)/2) * tempSigma - (nrow(x)/2) * tempDelta +
        0.5 * tempsum + lambda * tempmus + alpha * (ncol(x)/2) *
        tempSigma2 + beta * (nrow(x)/2) * tempDelta2
      return(obj)
    }

  if (is.null(Cs.init)) {
    Cs <- kmeans(x, k, nstart = 20)$cluster
  }
  else {
    Cs <- Cs.init
  }
  if (is.null(Ds.init)) {
    Ds <- kmeans(t(x), r, nstart = 20)$cluster
  }
  else {
    Ds <- Ds.init
  }
  if (center == TRUE) {
    mustemp <- mean(x)
    x <- x - mustemp
  }
  cl <- match.call()
  Cslist <- list()
  Dslist <- list()
  if (is.null(Sigma.init)) {
    Sigma <- diag(1, nrow = nrow(x), ncol = nrow(x))
    Delta <- diag(1, nrow = ncol(x), ncol = ncol(x))
    mus <- sparseBC:::updateMusMatrix(x, Cs, Ds, lambda, Sigma, Delta)
    objs <- 1e+15
    improvement <- 1e+10
    i <- 1
    while (improvement > (threshold) && i <= max.iter) {
      cat(i)
      mus <- sparseBC:::updateMusMatrix(x, Cs, Ds, lambda, Sigma,
                             Delta)
      Sigma <- sparseBC:::updateSigma(Delta, mus[Cs, Ds], x, alpha,
                           Cs, Ds)
      Delta <- sparseBC:::updateDelta(Sigma, mus[Cs, Ds], x, beta,
                           Cs, Ds)
      Cs <- sparseBC:::UpdateRowCluster(x, Sigma, Delta, mus[Cs, Ds],
                             Cs, Ds)
      Cs <- sparseBC:::ReNumberMatrix(Cs)
      Cslist[[i]] <- Cs
      mus <- sparseBC:::updateMusMatrix(x, Cs, Ds, lambda, Sigma,
                             Delta)
      Sigma <- sparseBC:::updateSigma(Delta, mus[Cs, Ds], x, alpha,
                           Cs, Ds)
      Delta <- sparseBC:::updateDelta(Sigma, mus[Cs, Ds], x, beta,
                           Cs, Ds)
      Ds <- sparseBC:::UpdateColumnCluster(x, Sigma, Delta, mus[Cs,
                                                     Ds], Cs, Ds)
      Ds <- sparseBC:::ReNumberMatrix(Ds)
      objs <- c(objs, MatrixObjective2(x, mus[Cs, Ds], Cs,
                                      Ds, Sigma, Delta, lambda, alpha, beta))
      Dslist[[i]] <- Ds
      i <- i + 1
      improvement <- abs(objs[i] - objs[i - 1])
    }
  }
  if (min(diff(objs)) > 0)
    print("Warning: objective values decreases")
  if (center == TRUE) {
    mus <- mus + mustemp
  }
  out <- list()
  class(out) <- "matrixBC"
  out$Cs <- Cslist[[i - 1]]
  out$Ds <- Dslist[[i - 1]]
  out$objs <- objs
  out$mus <- mus[out$Cs, out$Ds]
  out$Mus <- mus
  out$Sigma <- Sigma
  out$Delta <- Delta
  out$iteration <- i
  out$cl <- cl
  return(out)
}
