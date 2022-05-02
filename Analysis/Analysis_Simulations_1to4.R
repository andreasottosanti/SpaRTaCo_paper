#--------------------------------------------------
# Analyse Simulation Results (Simulations 1,...,4)
#--------------------------------------------------

library(spartaco)
library(sparseBC)
library(blockcluster)
# The file contains the MVNb model with the slight modification we implementend.
# See Section 4.2 for more details
source(".../matrixBC_modified.R")

Scenario <- 1                  # simulation number
data.directory <- "..."        # substitute '...' with the path where the data are stored
results.directory <- "..."     # substitute '...' with the path where the results are stored
n_replicas <- 10               # the number of replicas of the experiment

loglikelihoods <- rep(-Inf, n_replicas)

# ---Create the matrices where the clustering labels are stored
Best.Cs <- matrix(0,600,n_replicas)
Best.Cs.sparseBC <- array(0, c(600, n_replicas, 3))
Best.Cs.BC <- matrix(0,600,n_replicas)
Best.Cs.kmeans <- matrix(0,600,n_replicas)
Best.Cs.LBM <- matrix(0,600,n_replicas)
Best.Cs.MVNb <- array(0, c(600, n_replicas, 3))
Best.Ds <- matrix(0,600,n_replicas)
Best.Ds.sparseBC <- array(0, c(600, n_replicas, 3))
Best.Ds.BC <- matrix(0,600,n_replicas)
Best.Ds.kmeans <- matrix(0,600,n_replicas)
Best.Ds.LBM <- matrix(0,600,n_replicas)
Best.Ds.MVNb <- array(0, c(600, n_replicas, 3))

# ---Create the matrices where the CER values are stored
CER.values <- matrix(0,n_replicas,2)
CER.sparseBC <- array(0, c(n_replicas, 2, 3))
CER.BC <- matrix(0,n_replicas,2)
CER.kmeans <- matrix(0,n_replicas,2)
CER.LBM <- matrix(0,n_replicas,2)
CER.MVNb <- array(0, c(n_replicas, 2, 3))

for(i in 1:n_replicas){
  # ---load the simulation replica
  load(paste(data.directory,"/Scenario",Scenario,"_",i,".Rdata",sep=""))

  # ---load the parallel runs of SpaRTaCo on the same replica
  file.list <- dir(results.directory)[grepl(paste("Scenario",Scenario,"_",i,"_",sep=""), dir(results.directory))]
  best.spartaco <- CombineSpartaco(x = file.list, search.dir = results.directory, compute.uncertainty = F)
  CER.values[i,1] <- CER(Simulation$original.Cs, best.spartaco$Cs)
  CER.values[i,2] <- CER(Simulation$original.Ds, best.spartaco$Ds)

  # ---sparseBC
  sparseBC1 <- sparseBC(Simulation$x, k = model.dim[1], r = model.dim[2], lambda = 1)
  sparseBC2 <- sparseBC(Simulation$x, k = model.dim[1], r = model.dim[2], lambda = 10)
  sparseBC3 <- sparseBC(Simulation$x, k = model.dim[1], r = model.dim[2], lambda = 20)
  Best.Cs.sparseBC[,i,1] <- sparseBC1$Cs
  Best.Ds.sparseBC[,i,1] <- sparseBC1$Ds
  Best.Cs.sparseBC[,i,2] <- sparseBC2$Cs
  Best.Ds.sparseBC[,i,2] <- sparseBC2$Ds
  Best.Cs.sparseBC[,i,3] <- sparseBC3$Cs
  Best.Ds.sparseBC[,i,3] <- sparseBC3$Ds
  CER.sparseBC[i,1,1] <- CER(Simulation$original.Cs, sparseBC1$Cs)
  CER.sparseBC[i,2,1] <- CER(Simulation$original.Ds, sparseBC1$Ds)
  CER.sparseBC[i,1,2] <- CER(Simulation$original.Cs, sparseBC2$Cs)
  CER.sparseBC[i,2,2] <- CER(Simulation$original.Ds, sparseBC2$Ds)
  CER.sparseBC[i,1,3] <- CER(Simulation$original.Cs, sparseBC3$Cs)
  CER.sparseBC[i,2,3] <- CER(Simulation$original.Ds, sparseBC3$Ds)

  # ---BC
  BC <- sparseBC(Simulation$x, k = model.dim[1], r = model.dim[2], lambda = 0)
  CER.BC[i,1] <- CER(Simulation$original.Cs, BC$Cs)
  CER.BC[i,2] <- CER(Simulation$original.Ds, BC$Ds)
  Best.Cs.BC[,i] <- BC$Cs
  Best.Ds.BC[,i] <- BC$Ds

  # ---kmeans
  kmrow <- kmeans(Simulation$x, centers = model.dim[1], nstart = 10)
  kmcol <- kmeans(t(Simulation$x), centers = model.dim[2], nstart = 10)
  CER.kmeans[i,1] <- CER(Simulation$original.Cs, kmrow$clust)
  CER.kmeans[i,2] <- CER(Simulation$original.Ds, kmcol$clust)
  Best.Cs.kmeans[,i] <- kmrow$clust
  Best.Ds.kmeans[,i] <- kmcol$clust

  # ---LBM
  lbm <- coclusterContinuous(data = Simulation$x, nbcocluster = model.dim)
  CER.LBM[i,1] <- CER(Simulation$original.Cs, lbm@rowclass)
  CER.LBM[i,2] <- CER(Simulation$original.Ds, lbm@colclass)
  Best.Cs.LBM[,i] <- lbm@rowclass
  Best.Ds.LBM[,i] <- lbm@colclass

  # ---MVNb
  MVNb1 <- matrixBC2(x = Simulation$x, k = model.dim[1], r = model.dim[2], lambda = 1, alpha = .25, beta = .25)
  MVNb2 <- matrixBC2(x = Simulation$x, k = model.dim[1], r = model.dim[2], lambda = 10, alpha = 2.5, beta = 2.5)
  MVNb3 <- matrixBC2(x = Simulation$x, k = model.dim[1], r = model.dim[2], lambda = 20, alpha = 5, beta = 5)
  CER.MVNb[i,1,1] <- CER(Simulation$original.Cs, MVNb1$Cs)
  CER.MVNb[i,2,1] <- CER(Simulation$original.Ds, MVNb1$Ds)
  CER.MVNb[i,1,2] <- CER(Simulation$original.Cs, MVNb2$Cs)
  CER.MVNb[i,2,2] <- CER(Simulation$original.Ds, MVNb2$Ds)
  CER.MVNb[i,1,3] <- CER(Simulation$original.Cs, MVNb3$Cs)
  CER.MVNb[i,2,3] <- CER(Simulation$original.Ds, MVNb3$Ds)
  Best.Cs.MVNb[,i,1] <- MVNb1$Cs
  Best.Ds.MVNb[,i,1] <- MVNb1$Ds
  Best.Cs.MVNb[,i,2] <- MVNb2$Cs
  Best.Ds.MVNb[,i,2] <- MVNb2$Ds
  Best.Cs.MVNb[,i,3] <- MVNb3$Cs
  Best.Ds.MVNb[,i,3] <- MVNb3$Ds
}

