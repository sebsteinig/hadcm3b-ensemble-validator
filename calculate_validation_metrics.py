import os
import logging
import argparse

from metrics import global_productivity_fluxes

# user settings
username = os.getlogin()
# storage on BRIDGE servers
data_dir = f"/export/silurian/array-01/{username}/ensembles"

# set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('calculate_validation_metrics.log'),
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
    
    for id in ids:
        output_dir = os.path.join(data_dir, id, "processed")
        create_directory(output_dir)
        
        for metric in metrics:
            logging.info(f"Calculating metric: {metric}")
            metric_dir = os.path.join(output_dir, metric)
            create_directory(metric_dir)

            if metric == "global_productivity_fluxes":
                global_productivity_fluxes(data_dir, id, output_dir, logging)

            
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Calculate validation metrics from UM ensemble output files.')
    parser.add_argument('experiments', nargs='+', help='List of experiment names to process')
    
    args = parser.parse_args()
    
    for experiment in args.experiments:
        logging.info(f"Metric calculation script started for experiment: {experiment}")
        main(experiment)
        logging.info(f"Metric calculation script finished for experiment: {experiment}")
