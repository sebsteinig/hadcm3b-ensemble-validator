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
clim_end_year = 1871

# set up logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[
        logging.FileHandler("plot_validation_metrics.log"),
        logging.StreamHandler(),
    ],
)

# metrics = ["global_productivity_fluxes", "global_carbon_stores", "RECCAP_stores_vs_fluxes"]
# metrics = ["global_productivity_fluxes", "global_carbon_stores", "RECCAP_stores_vs_fluxes", "overview_table", "PFT_maps"]
# metrics = ["global_productivity_fluxes", "global_carbon_stores", "RECCAP_stores_vs_fluxes"]
# metrics = ["global_veg_fractions"]
# metrics = ["overview_table", "RECCAP_stores_vs_fluxes"]
metrics = ["overview_table"]
# metrics = ["PFT_maps"]
# metrics = ["RECCAP_stores_vs_fluxes"]


# these ensemble members will be highlighted in the scatter plots
# this is also the list of IDs used in the PFT maps
# highlight_ids = {"xqaba": "control", "xqabs": "#1 score", "xqabV": "#2 score", "xQabu": "#3 score", "XqabR": "#4 score", "Xqabb": "#5 score", "XqAbw": "#6 score", "xQabn": "#7 score", "xqabQ": "#8 score", "XqAbf": "#9 score"}
# highlight_ids = {"xqaqa": "control"}
# highlight_ids = {"xqacz": "#1 score", "Xqacs": "#2 score", "Xqacp": "#3 score", "xqaca" : "control"}
# highlight_ids = {"xqaqa": "control", "XqAqd": "#1 score", "xqaqq": "#2 score", "xQaqc": "#3 score", "xqaQs": "#4 score", "Xqaqp": "#5 score", "XqaqW": "#6 score", "xqaqg": "#7 score", "xqaqC": "#8 score", "xqaqv": "#9 score"}
# highlight_ids = {"xqaua": "control"}
# highlight_ids = {"xqara": "control", "Xqarh": "#1 score", "Xqarp": "#2 score", "xqare": "#3 score", "xqaRw": "#4 score", "XqaRq": "#5 score", "xqarR": "#6 score", "xqarp": "#7 score", "Xqaru": "#8 score", "XQarl": "#9 score"}
# highlight_ids = {"xqaua": "control", "Xqaul": "#1 score", "xqAuw": "#2 score", "Xqaum": "#3 score", "xqauO": "#4 score", "xqaUu": "#5 score", "xqaub": "#6 score", "xQaup": "#7 score", "xqauj": "#8 score", "xqaUs": "#9 score"}
# highlight_ids = {"XqAqd": "#1 score", "xqaqa": "control", "XQari": "#256 score"}
# extended
# highlight_ids = {"xqaqa": "control", "Xqarh": "#1 score", "XqAqs": "#2 score", "Xqaul": "#3 score", "Xqarp": "#4 score", "xqare": "#5 score", "XqAuf": "#6 score", "xqaQs": "#7 score", "xqaRw": "#8 score", "xQaqc": "#9 score",
#                  "Xqaul": "#10 score", "Xqaqp": "#11 score", "XqaRq": "#12 score", "xQaup": "#13 score", "xqAuw": "#14 score", "XqAqg": "#15 score", "xqarR": "#16 score", "XqaqW": "#17 score", "XqAux": "#18 score", "XqAuc": "#19 score", "xqaqg": "#20 score",
#                  "xqaQh": "low VCRIT #1", "xqaRd": "low VCRIT #2", "XqArn": "low VCRIT #3"}
# reduced
highlight_ids = {"xqaqa": "control", "XqaqW": "#1 score", "Xqarh": "#2 score", "XqAuc": "#3 score", "xQaup": "#4 score", "XqaRq": "#5 score", "XqAuf": "#6 score", "Xqaul": "#7 score", "xqare": "#8 score", "xQaqc": "#9 score",
                 "XqAqs": "#10 score", "Xqarp": "#11 score", "Xqaqp": "#12 score", "xqaRw": "#13 score", "xqarE": "#14 score", "XQarc": "#15 score", "XqAqg": "#16 score", "XqauI": "#17 score", "xqaqg": "#18 score", "xqaQs": "#19 score", "xqauj": "#20 score",
                 "xqaQh": "low VCRIT #1", "xqaRd": "low VCRIT #2", "XqArn": "low VCRIT #3"}
# metrics to include in overview table and skill score calculation
selected_metrics = {
    "GPP": {"weight": 8},
    "CVeg": {"weight": 8},
    # "CSoil": {"weight": 1},
    "Tr30SN": {"weight": 4},
    "Tr30-90N": {"weight": 4},
    "AMZTrees": {"weight": 4},
    "GM_BL": {"weight": 2},
    "GM_NL": {"weight": 2},
    "GM_C3": {"weight": 2},
    "GM_C4": {"weight": 2},
    "GM_BS": {"weight": 2},
    "rmse_BL": {"weight": 2},
    "rmse_NL": {"weight": 2},
    "rmse_C3": {"weight": 2},
    "rmse_C4": {"weight": 2},
    "rmse_BS": {"weight": 2},
    # "RECCAP_VEG_C": {"weight": 1},
    # "RECCAP_GPP": {"weight": 1},
    # "RECCAP_SOIL_C": {"weight": 1},
    # "RECCAP_TAU": {"weight": 1},
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
