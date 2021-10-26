library(spartaco)

load(paste(data.directory,".../LIBD_subj151673_1585spots.Rdata",sep=""))  # substitute '...' with the path where the data are stored

results.directory <- "..."     # substitute '...' with the path where the results will be stored

K.values <- 1:3  # range of row clusters
R.values <- 4:8  # range of column clusters
n_starting_points <- 5  # This value gives the number of times you want to run the estimation algorithm in parallel

parallel::mclapply(K.values, function(k) {
  parallel::mclapply(R.values, function(r){
    parallel::mclapply(1:n_starting_points, function(j){
      results <- spartaco(x = Data$x,
                      coordinates = Data$coordinates,
                      K = k, R = r,
                      traceRatio = 10, max.iter = 5*10^3,
                      metropolis.iterations = 150,
                      estimate.iterations = 10)
      save(results, file = paste(results.directory,"/spatialLIBD_K",k,"_R",r,"_run",j,"_",
                                 substr(Sys.time(),1,10),"_",
                                 substr(Sys.time(),12,13),"_",
                                 substr(Sys.time(),15,16),"_",
                                 substr(Sys.time(),18,19),".Rdata",sep=""))
      }, mc.cores = 1)
    }, mc.cores = 1)
}, mc.cores = 1)

# You can set the mc.cores parameters according to the characteristics of your machine
