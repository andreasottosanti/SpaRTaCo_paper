library(spartaco)
library(microbenchmark)

Scenario <- 1
data.directory <- "..."        # substitute '...' with the path where the data are stored
results.directory <- "..."     # substitute '...' with the path where the results will be stored

n_replicas <- 10  # This depends from how many replicas of the experiment you generated
n_starting_points <- 5  # This value gives the number of times you want to run the estimation algorithm in parallel

parallel::mclapply(1:n_replicas, function(i) {
  parallel::mclapply(1:n_starting_points, function(j){
    print(paste("Loading Scenario",Scenario,"_",i,sep=""))
    load(paste(data.directory,"/Scenario",Scenario,"_",i,".Rdata",sep=""))
    file.name = paste(results.directory,"/Scenario",Scenario,"_K3_R3_run",j,".Rdata",sep="")
    exec.time <- microbenchmark(
      results <- spartaco(x = x, coordinates = coordinates, K = 3, R = 3,
                          max.iter = 3*10^3,
                          metropolis.iterations = 150,
                          estimate.iterations = 100, verbose=T,
                          conv.criterion = list(iterations = 30, epsilon = 1e-04),
                          save.options = list(after = 100, file.name = file.name)),
      times = 1)
    results$exec.time <- exec.time$time
    save(results, file = file.name)
  },
  mc.cores = 5)
}, mc.cores = 5)

# You can set the mc.cores parameters according to the characteristics of your machine
# Here we set 25 cores
