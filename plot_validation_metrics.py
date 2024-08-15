import os
import logging
import argparse

from common import parse_jobs, create_directory, get_table_metrics
from plotting import (
    plot_timeseries,
    plot_parameter_scatter,
    plot_RECCAP_stores_vs_fluxes,
    plot_overview_table,
    plot_skill_score_scatter,
    plot_PFT_maps,
)

# user settings
username = os.getlogin()
# storage on BRIDGE servers
data_dir = f"/export/silurian/array-01/{username}/ensembles"

clim_start_year = 1859
clim_end_year = 1890

# set up logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[
        logging.FileHandler("plot_validation_metrics.log"),
        logging.StreamHandler(),
    ],
)

# metrics = ["RECCAP_stores_vs_fluxes"]
# metrics = ["global_productivity_fluxes", "global_carbon_stores", "RECCAP_stores_vs_fluxes", "overview_table", "PFT_maps"]
# metrics = ["RECCAP_stores_vs_fluxes", "overview_table", "PFT_maps"]
# metrics = ["global_veg_fractions"]
# metrics = ["RECCAP_stores_vs_fluxes"]
metrics = ["overview_table"]
# metrics = ["PFT_maps"]

# these ensemble members will be highlighted in the scatter plots
# this is also the list of IDs used in the PFT maps
highlight_ids = {"xQabn": "#1 score", "xQaby": "#2 score", "XqAbw": "#4 score", "xpwca": "acang", "xpwcn" : "no mods control"}
# highlight_ids = {"xpzna": "control"}

# metrics to include in overview table and skill score calculation
selected_metrics = {
    "GPP": {"weight": 5},
    "CVeg": {"weight": 5},
    "CSoil": {"weight": 5},
    "Tr30SN": {"weight": 3},
    "Tr30-90N": {"weight": 3},
    "BareSoil": {"weight": 3},
    "RECCAP_VEG_C": {"weight": 1},
    "RECCAP_GPP": {"weight": 1},
    "RECCAP_SOIL_C": {"weight": 1},
    "RECCAP_TAU": {"weight": 1},
}
table_metrics = get_table_metrics(selected_metrics)


def main(experiment):
    model_params = parse_jobs(f"./id_lists/{experiment}_parameters.json")
    logging.info(f"Read {len(model_params)} IDs from log file")

    output_dir = f"./ensemble_validation_plots"
    create_directory(output_dir, logging)

    # To Do: create combined job list from all experiments

    # create one output file for each metric
    for metric in metrics:
        if metric in [
            "global_productivity_fluxes",
            "global_carbon_stores",
            "global_veg_fractions",
        ]:
            plot_timeseries(
                metric, model_params, data_dir, experiment, output_dir, logging
            )
            plot_parameter_scatter(
                metric,
                model_params,
                data_dir,
                experiment,
                output_dir,
                logging,
                clim_start_year,
                clim_end_year,
                highlight_ids
            )
        elif metric == "RECCAP_stores_vs_fluxes":
            plot_RECCAP_stores_vs_fluxes(
                model_params,
                data_dir,
                experiment,
                output_dir,
                logging,
                clim_start_year,
                clim_end_year,
                highlight_ids
            )
        elif metric == "overview_table":
            plot_overview_table(
                model_params,
                data_dir,
                experiment,
                output_dir,
                logging,
                table_metrics,
                clim_start_year,
                clim_end_year,
            )
            plot_skill_score_scatter(
                model_params,
                experiment,
                output_dir,
                logging,
                clim_start_year,
                clim_end_year,
                highlight_ids
            )
        elif metric == "PFT_maps":
            plot_PFT_maps(
                model_params,
                data_dir,
                experiment,
                output_dir,
                logging,
                clim_start_year,
                clim_end_year,
                highlight_ids
            )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Plot validation metrics from UM ensemble output files."
    )
    parser.add_argument(
        "experiments", nargs="+", help="List of experiment names to include into plots"
    )

    args = parser.parse_args()

    logging.info(
        f"Metric plotting script started for experiment(s): {args.experiments[0]}"
    )
    main(args.experiments[0])
    logging.info(
        f"Metric plotting script finished for experiment(s): {args.experiments[0]}"
    )
