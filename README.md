# Simulations and Data Analysis using SpaRTaCo
This repository contains the scripts to reproduce the simulations and the analysis shown in the paper '*Co-clustering of Spatially Resolved Transcriptomic Data*' by Andrea Sottosanti and Davide Risso available [here](https://arxiv.org/abs/2110.04872).

To run SpaRTaCo, you need to install the library `spartaco` available [here](https://github.com/andreasottosanti/spartaco).

&nbsp;
## Reproduce the simulation experiments
To reproduce the spatial experiments proposed in the article, run the files contained into the directory `.../Create_Simulations`. Every file generates and saves multiple replicas of the same experiment. You can choose the number of replicas with `n_replicas` (default is `n_replicas = 10` for Simulations 1,..,4 and `n_replicas = 1` for Simulation 5). You need also to specify a directory where to save the simulated datasets.

&nbsp;
## Run SpaRTaCo on the simulation experiments
SpaRTaCo can be run multiple times in parallel on each and every replica of the same experiment with just one script, `.../Analysis/ModelFitting_Simulations.R`.
Before running the code, 

1. set the number of the simulation scenario `Scenario`;
2. set `data.directory` with the path of the directory where you stored the replicas of the experiment;
3. set `results.directory` with the path of the directory where you the results to be saved;
4. set `n_replicas` equal to the number of replicas of the experiment you have generated;
5. set how many times you want to run the estimation algorithm in parallel with `n_starting_points`;
6. set the model options in the function `spartaco` (see the help for more details);
7. set the number of cores according to your machine characteristics.
