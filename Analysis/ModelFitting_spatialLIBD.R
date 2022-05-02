library(spartaco)
library(microbenchmark)

load(paste(data.directory,".../LIBD_subj151673.Rdata",sep=""))  # substitute '...' with the path where the data are stored

results.directory <- "..."     # substitute '...' with the path where the results will be stored

K.values <- 2  # range of row clusters
R.values <- 7:12  # range of column clusters
n_starting_points <- 5  # This value gives the number of times you want to run the estimation algorithm in parallel

parallel::mclapply(K.values, function(k) {
  parallel::mclapply(R.values, function(r){
    parallel::mclapply(1:n_starting_points, function(j){
      file.name = paste(results.directory,"/LIBD151673_K",k,"_R",r,"_run",j,".Rdata",sep="")
      exec.time <- microbenchmark(
        results <- spartaco(x = x, coordinates = coordinates, K = k, R = r,
                            max.iter = 3*10^3,
                            metropolis.iterations = 150,
                            estimate.iterations = 100, verbose=T,
                            conv.criterion = list(iterations = 30, epsilon = 1e-04),
                            save.options = list(after = 100, file.name = file.name)),
        times = 1)
      results$exec.time <- exec.time$time
      save(results, file = file.name)
      }, mc.cores = 1)
    }, mc.cores = 1)
}, mc.cores = 1)

# You can set the mc.cores parameters according to the characteristics of your machine
