import os
import logging
import argparse

from common import parse_jobs, create_directory
from plotting import plot_timeseries, plot_parameter_scatter

# user settings
username = os.getlogin()
# storage on BRIDGE servers
data_dir = f"/export/silurian/array-01/{username}/ensembles"

clim_start_year = 1860
clim_end_year   = 1870

# set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('plot_validation_metrics.log'),
        logging.StreamHandler()
    ]
)

# metrics = ["global_productivity_fluxes"]
# metrics = ["global_productivity_fluxes", "global_carbon_stores"]
metrics = ["global_carbon_stores"]

    
def main(experiment):
    model_params = parse_jobs(f"./id_lists/{experiment}_parameters.json")
    logging.info(f"Read {len(model_params)} IDs from log file")
    
    output_dir = f"./ensemble_validation_plots"
    create_directory(output_dir, logging)

    # To Do: create combined job list from all experiments
    
    # create one output file for each metric
    for metric in metrics:
        if metric in [ "global_productivity_fluxes" ]:
            plot_timeseries(metric, model_params, data_dir, experiment, output_dir, logging)
            plot_parameter_scatter(metric, model_params, data_dir, experiment, output_dir, logging, clim_start_year, clim_end_year)
        if metric in [ "global_carbon_stores" ]:
            plot_timeseries(metric, model_params, data_dir, experiment, output_dir, logging)
            plot_parameter_scatter(metric, model_params, data_dir, experiment, output_dir, logging, clim_start_year, clim_end_year)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Plot validation metrics from UM ensemble output files.')
    parser.add_argument('experiments', nargs='+', help='List of experiment names to include into plots')
    
    args = parser.parse_args()
    
    logging.info(f"Metric plotting script started for experiment(s): {args.experiments[0]}")
    main(args.experiments[0])
    logging.info(f"Metric plotting script finished for experiment(s): {args.experiments[0]}")
