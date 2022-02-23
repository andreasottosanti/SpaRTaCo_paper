#----------------------------------------
### BUILD THE SIMULATION EXPERIMENT 0 ###
#----------------------------------------

# ---import the coordinates a create the spot labels
load("~/Create_Simulations/coordinates.Rdata")
Ds <- c(rep(1,200),rep(2,200),rep(3,200))  # this is the vector of the column clustering labels
tableDs <- as.vector(table(Ds))
plot(coordinates, col = Ds, pch = 16)
Dist <- as.matrix(stats::dist(coordinates))

# ---prepare the vector of row clustering labels
n1 <- n2 <- n3 <- 200  # size of the row clusters
Cs <- numeric(n1 + n2 + n3)   # this is the vector of the row clustering labels
Cs[1:n1] <- 1
Cs[(n1+1):(n1+n2)] <- 2
Cs[(n1+n2+1):(n1+n2+n3)] <- 3

# ---define the spatial signal-to-noise ratios
ratio <- 3
Ratios <- matrix(ratio,3,3)
# convert the spatial signal-to-noise ratios into the values of tau and xi
c_Delta <- 10    # this is the quantity tau[k,r]+xi[k,r]
Tau <- c_Delta * Ratios/(Ratios+1)
Xi <- c_Delta - Tau

# ---define the spatial covariance parameters
KernelParameters <- list(Exponential = 50,
                         RatQuadr = c(50, 2),
                         Gaussian = 70)

n_replicas <- 10   # set the number of replicas (10 in the paper)
for(j in 1:n_replicas){
  # create the empty data matrix
  x <- matrix(0, sum(n1+n2+n3), nrow(coordinates))

  # ---define the mean parameters
  mu <- 5
  Mu <- matrix(c(0,mu,-mu,-mu,0,mu,mu,-mu,0),3,3)

  # ---draw the Sigma2 matrices
  Sigma <- array(0, c(n1, n1, 3))
  Sigma[,,1] <- rWishart(n = 1, df = n1+10, Sigma = diag(.03, n1))[,,1]
  Sigma[,,3] <- rWishart(n = 1, df = n1, Sigma = Sigma[,,1]/150)[,,1]
  Sigma[,,2] <- rWishart(n = 1, df = n1, Sigma = Sigma[,,1]/150)[,,1]
  Sigma[,,1] <- rWishart(n = 1, df = n1, Sigma = Sigma[,,1]/150)[,,1]

  # ---draw the blocks of x
  Kernels <- array(0, c(sum(Ds == 1), sum(Ds == 1), 3))
  for(r in 1:3){
    Kernels[,,r] <- exp(-Dist[Ds == r, Ds == r]^2/(2*KernelParameters$Gaussian^2))
    for(k in 1:3){
      Delta <- Tau[k,r]*Kernels[,,r]+diag(Xi[k,r], nrow = tableDs[r])
      x[Cs == k, Ds == r] <- Mu[k, r] + t(chol(Sigma[,,k])) %*% matrix(rnorm(sum(Cs == k)*tableDs[r]), sum(Cs == k), tableDs[r]) %*% chol(Delta)
    }
  }

  # ---save the data
  Simulation <- list(x = x,
                     Mu = Mu,
                     original.Cs = Cs,
                     original.Ds = Ds,
                     coordinates = coordinates,
                     Tau = Tau,
                     Xi = Xi,
                     KernelParameters = KernelParameters,
                     Sigma = Sigma)
  save(Simulation, file = paste("~/Scenario0_",j,".Rdata",sep=""))
}
