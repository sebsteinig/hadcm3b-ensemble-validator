import os
import logging
import argparse

from common import parse_jobs, create_directory
from metrics import global_productivity_fluxes, global_carbon_stores, global_veg_fractions_model, global_veg_fractions_obs

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

# metrics = ["global_veg_fractions"]
metrics = ["global_productivity_fluxes", "global_carbon_stores", "global_veg_fractions"]
# metrics = ["global_carbon_stores"]
# metrics = ["global_productivity_fluxes"]
# metrics = ["global_veg_fractions"]



def main(experiment):
    model_params = parse_jobs(f"./id_lists/{experiment}_parameters.json")
    logging.info(f"Read {len(model_params)} IDs from log file")
    
    # process each ensemble member
    for id, params in model_params.items():
        output_dir = os.path.join(data_dir, id, "processed")
        create_directory(output_dir, logging)
        
        for metric in metrics:
            logging.info(f"Calculating metric: {metric}")
            metric_dir = os.path.join(output_dir, metric)
            create_directory(metric_dir, logging)

            if metric == "global_productivity_fluxes":
                global_productivity_fluxes(data_dir, id, output_dir, logging)
            elif metric == "global_carbon_stores":
                global_carbon_stores(data_dir, id, output_dir, logging)
            elif metric == "global_veg_fractions":
                global_veg_fractions_model(data_dir, id, output_dir, logging)

    # process some reference observations
    for metric in metrics:
        if metric == "global_veg_fractions":
            global_veg_fractions_obs('./observations/qrparm.veg.frac_igbp.pp.nc', logging)

            

            
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Calculate validation metrics from UM ensemble output files.')
    parser.add_argument('experiments', nargs='+', help='List of experiment names to process')
    
    args = parser.parse_args()
    
    for experiment in args.experiments:
        logging.info(f"Metric calculation script started for experiment: {experiment}")
        main(experiment)
        logging.info(f"Metric calculation script finished for experiment: {experiment}")
