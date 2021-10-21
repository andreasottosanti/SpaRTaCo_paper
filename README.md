# Simulations and Data Analysis using SpaRTaCo
This repository contains the scripts to reproduce the simulations and the analysis shown in the paper '*Co-clustering of Spatially Resolved Transcriptomic Data*' by Andrea Sottosanti and Davide Risso available [here](https://arxiv.org/abs/2110.04872).

To run SpaRTaCo, you need to install the library `spartaco` available [here](https://github.com/andreasottosanti/spartaco).


## Reproduce the simulation experiments
To reproduce the spatial experiments proposed in the article, you can run the files contained into the directory `.../Create_Simulations`. Every file will generate and save multiple replicas of the same experiment. You can choose the number of replicas with `n_replicas` (default is `n_replicas = 10` for Simulations 1,..,4 and `n_replicas = 1` for Simulation 5). You need also to specify a directory where you can save the simulated datasets.
