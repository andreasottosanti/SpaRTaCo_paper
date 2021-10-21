#----------------------------------------
### BUILD THE SIMULATION EXPERIMENT 5 ###
#----------------------------------------

# ---import the coordinates a create the spot labels
load("~/Create_Simulations/coordinates.Rdata")
Ds <- c(rep(1,200),rep(2,200),rep(3,200))  # this is the vector of the column clustering labels
tableDs <- as.vector(table(Ds))
plot(coordinates, col = Ds, pch = 16)
Dist <- as.matrix(stats::dist(coordinates))


# ---prepare the vector of row clustering labels (see Figure 4 (d) of the manuscript)
N <- matrix(c(rep(200,6),100,200,300),3,3)
Cs <- matrix(0, 600, 3)
Cs[1:200,] <- 1
Cs[(201):(400),] <- 2
Cs[(401):(600),] <- 3
Cs[101:200,2] <- Cs[301:400,2] <- 2
Cs[201:300,2] <- Cs[501:600,2] <- 3
Cs[401:500,2] <- 1
Cs[101:200,3] <- 3
# this is the new row clustering labelling system introduced to evaluate the results
subCluster.Cs <- numeric(nrow(Cs))
subCluster.Cs[Cs[,1] == 1 & Cs[,2] == 1 & Cs[,3] == 1] <- 1
subCluster.Cs[Cs[,1] == 1 & Cs[,2] == 2 & Cs[,3] == 3] <- 2
subCluster.Cs[Cs[,1] == 2 & Cs[,2] == 3 & Cs[,3] == 2] <- 3
subCluster.Cs[Cs[,1] == 2 & Cs[,2] == 2 & Cs[,3] == 2] <- 4
subCluster.Cs[Cs[,1] == 3 & Cs[,2] == 1 & Cs[,3] == 3] <- 5
subCluster.Cs[Cs[,1] == 3 & Cs[,2] == 3 & Cs[,3] == 3] <- 6



# ---define the spatial signal-to-noise ratios 
Ratios <- matrix(0,3,3)
Ratios[1,] <- c(0, 3, 1)
Ratios[2,] <- c(1, 0, 3)
Ratios[3,] <- c(3, 1, 0)
# convert the spatial signal-to-noise ratios into the values of tau and xi
c_Delta <- 10    # this is the quantity tau[k,r]+xi[k,r]
Tau <- c_Delta * Ratios/(Ratios+1)
Xi <- c_Delta - Tau

# ---define the spatial covariance parameters
KernelParameters <- list(Exponential = 50,
                         RatQuadr = c(50, 2),
                         Gaussian = 70)

n_replicas <- 1   # set the number of replicas (1 in the paper)
for(j in 1:n_replicas){
  # create the empty data matrix
  x <- matrix(0, 600, nrow(coordinates))
  
  # ---draw the Sigma2 matrices
  Sigma <- list()
  for(r in 1:3){
    Sigma[[r]] <- matrix(0, nrow(x), nrow(x))
    Sigma[[r]][Cs[,r] == 1,Cs[,r] == 1] <- rWishart(n = 1, df = N[1,r]+10, Sigma = diag(.03, N[1,r]))[,,1]
    Sigma[[r]][Cs[,r] == 2,Cs[,r] == 2] <- rWishart(n = 1, df = N[2,r]+30, Sigma = diag(.05, N[2,r]))[,,1]
    Sigma[[r]][Cs[,r] == 3,Cs[,r] == 3] <- rWishart(n = 1, df = N[3,r], Sigma = rWishart(n = 1, df = N[3,r]+10, Sigma = diag(.03, N[3,r]))[,,1]/150)[,,1]
  }
  
  # ---draw the blocks of x
  Kernels <- array(0, c(sum(Ds == 1), sum(Ds == 1), 3))
  for(r in 1:3){
    if(r == 1) Kernels[,,r] <- exp(-Dist[Ds == r, Ds == r]/KernelParameters$Exponential)
    if(r == 2) Kernels[,,r] <- (1+Dist[Ds == r, Ds == r]^2/(2*KernelParameters$RatQuadr[2]*KernelParameters$RatQuadr[1]^2))^(-KernelParameters$RatQuadr[2])
    if(r == 3) Kernels[,,r] <- exp(-Dist[Ds == r, Ds == r]^2/(2*KernelParameters$Gaussian^2))
    for(k in 1:3){
      Delta <- Tau[k,r]*Kernels[,,r]+diag(Xi[k,r], nrow = tableDs[r])
      x[Cs[,r] == k, Ds == r] <- t(chol(Sigma[[r]][Cs[,r]==k,Cs[,r]==k])) %*% matrix(rnorm(sum(Cs[,r] == k)*tableDs[r]), sum(Cs[,r] == k), tableDs[r]) %*% chol(Delta)
    }
  }
    
  Simulation <- list(x = x, 
                     original.Cs = Cs,
                     subCluster.Cs = subCluster.Cs,
                     original.Ds = Ds,
                     coordinates = coordinates,
                     Tau = Tau,
                     Xi = Xi,
                     KernelParameters = KernelParameters,
                     Sigma = Sigma)
  save(Simulation, file = paste("~/Scenario5_",j,".Rdata",sep=""))
}
