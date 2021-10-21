# Simulations and Data Analysis using SpaRTaCo
This repository contains the scripts to reproduce the simulations and the analysis shown in the paper '*Co-clustering of Spatially Resolved Transcriptomic Data*' by Andrea Sottosanti and Davide Risso available [here](https://arxiv.org/abs/2110.04872).

To run SpaRTaCo, you need to install the library `spartaco` available [here](https://github.com/andreasottosanti/spartaco).



## Reproduce the simulation experiments
To reproduce the spatial experiments proposed in the article, you can run the files contained into the directory `.../Create_Simulations`. Every file will generate and save multiple replicas of the same experiment. You can choose the number of replicas with `n_replicas` (default is `n_replicas = 10` for Simulations 1,..,4 and `n_replicas = 1` for Simulation 5). You need also to specify a directory where you can save the simulated datasets.



## Run SpaRTaCo on the simulation experiments
You can run SpaRTaCo multiple times in parallel on each and every replica of the same experiment with the script `.../Analysis/ModelFitting_Simulations.R`.
Before running the code, 

1. set the number of the simulation scenario `Scenario`;
2. set `n_replicas` equal to the number of replicas you made of this simulation experiment;
3. set the number of parallel run of the estimation algorithm with `n_starting_points`;
4. write the path of the directory which contains the replicas of the experiment;
5. set the model;
6. set the directory where you want to save the results (and eventually the name of the saved file);
7. set the number of cores.
