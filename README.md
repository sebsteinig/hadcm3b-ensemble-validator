# HadCM3B land carbon cycle ensemble validator
This repo contains scripts to automatically download, process and visualise HadCM3B UM ensemble jobs. These ensembles are intended to test a large range of model parameters for the interactive land carbon cylce component added by Chris Jones. The goal is to run the short, standard preindustrial simulation with hundreds (or thousands) of different tuning parameters to ultimately identify one (or several) promising parameter sets to use as a baseline configuration for the new HadCM3BL version with a fully coupled carbon cycle. 

Creation and running of the ensembles can be done via https://github.com/sebsteinig/hadcm3b-ensemble-generator, this repo aims to extract relevant terrestrial carbon cycle metrics from the output data and to compare them with observational targets to ultimately generate a simple skill score to rank the performance of individual ensemble members. It is intended to be copied and run in the user home directory on the BRIDGE servers.

## example workflow
Let's assume we want to analyse an ensemble of 250 simulations, each with random values of key land carbon cycle parameters to assess the model sensitivity against the parameter choice. For this, we would need to go through the following steps:

0. dependens on: ensemble simulation output 
- successfully created and ran the ensemble with https://github.com/sebsteinig/hadcm3b-ensemble-generator on BC4

1. download data
- transfer model output from standard BC4 location to BRIDGE Silurian server with `python download_files.py <EXPID>`
- the script can be configured to transfer only selected UM output streams (e.g. "pi" for the TRIFFID output)
- the script can be run iteratively and will only transfer/process new files

2. calculate validation metrics
- once the data is downloaded, we can calculate custom validation metrics for each file with `python calculate_validation_metrics.py <EXPID>`
- examples are global sums of vegetation and soil carbon, carbon fluxes or mean vegetation PFT fractions
- the list of metrics can be configured within the script (e.g. `metrics = ["global_productivity_fluxes", "global_carbon_stores", "global_veg_fractions"]`), with the actual definitions and code avaialble in the `metrics.py` module
- this modular approach should make the script easily extendable to new diagnostics in the future
- results are written to disk to have them readily available for plotting to speed up computation

3. plot validation metrics
- finally, we can plot the validation metrics with `python plot_validation_metrics.py <EXPID>`
- again, the list of metrics to polot can be set within the script (`metrics = ["global_productivity_fluxes", "global_carbon_stores", "RECCAP_stores_vs_fluxes"]`)
- the plots themselves can be defined (and expanded in the future) with the `plotting.py` module

