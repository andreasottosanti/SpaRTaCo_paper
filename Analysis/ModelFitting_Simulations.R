library(spartaco)

Scenario <- 1  
data.directory <- "..."        # substitute '...' with the path where the data are stored
results.directory <- "..."     # substitute '...' with the path where the results will be stored

n_replicas <- 10  # This depends from how many replicas of the experiment you generated
n_starting_points <- 5  # This value gives the number of times you want to run the estimation algorithm in parallel

parallel::mclapply(1:n_replicas, function(i) {
  parallel::mclapply(1:n_starting_points, function(j){
    print(paste("Loading Scenario",Scenario,"_",i,sep=""))
    load(paste(data.directory,"/Scenario",Scenario,"_",i,".Rdata",sep=""))  
    results <- spartaco(x = Simulation$x, 
                    coordinates = Simulation$coordinates, 
                    K = 3, R = 3,
                    traceRatio = 10, max.iter = 5*10^3,
                    metropolis.iterations = 150,
                    estimate.iterations = 10)
    save(results, file = paste(results.directory,"/Scenario",Scenario,"_replica",i,"_run",j,"_",
                               substr(Sys.time(),1,10),"_",
                               substr(Sys.time(),12,13),"_",
                               substr(Sys.time(),15,16),"_",
                               substr(Sys.time(),18,19),".Rdata",sep=""))
  },
  mc.cores = 5)
}, mc.cores = 5)    

# You can set the mc.cores parameters according to the characteristics of your machine
# Here we set 25 cores
