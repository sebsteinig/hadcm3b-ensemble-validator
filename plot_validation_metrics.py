import os
import logging
import argparse

from plotting import plot_timeseries, plot_parameter_scatter

# user settings
username = os.getlogin()
# storage on BRIDGE servers
data_dir = f"/export/silurian/array-01/{username}/ensembles"

# set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('plot_validation_metrics.log'),
        logging.StreamHandler()
    ]
)


# read job IDs
def read_ids_from_log(file_path):
    with open(file_path, 'r') as file:
        return [line.strip() for line in file.readlines()]


def create_directory(path):
    if not os.path.exists(path):
        os.makedirs(path)
        logging.info(f"Created directory: {path}")
    else:
        logging.info(f"Directory already exists: {path}")



metrics = ["global_productivity_fluxes"]

        
def main(experiment):
    id_list = f"./id_lists/{experiment}.log"
    ids = read_ids_from_log(id_list)
    logging.info(f"Read {len(ids)} IDs from log file")
    
    output_dir = f"./ensemble_validation_plots"
    create_directory(output_dir)

    # To Do: create combined job list from all experiments
    
    # create one output file for each metric
    for metric in metrics:
        if metric in [ "global_productivity_fluxes" ]:
            plot_timeseries(metric, data_dir, experiment, output_dir, logging)
            plot_parameter_scatter(metric, data_dir, experiment, output_dir, logging)


            
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Plot validation metrics from UM ensemble output files.')
    parser.add_argument('experiments', nargs='+', help='List of experiment names to include into plots')
    
    args = parser.parse_args()
    
    logging.info(f"Metric plotiing script started for experiment(s): {args.experiments[0]}")
    main(args.experiments[0])
    logging.info(f"Metric plotiing script finished for experiment(s): {args.experiments[0]}")
