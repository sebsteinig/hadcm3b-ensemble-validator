import os
import json
import xarray as xr
import csv


def parse_jobs(filename):
    """
    Parses a JSON file and returns a dictionary with ensemble_id as keys
    and their corresponding model parameters as values.

    Args:
        filename (str): The path to the JSON file.

    Returns:
        dict: A dictionary with ensemble_id as keys and parameter entries as values.
    """
    with open(filename, "r") as file:
        data = json.load(file)

    result = {}
    for entry in data:
        ensemble_id = entry.get("ensemble_id")
        if ensemble_id:
            result[ensemble_id] = entry

    return result


def create_directory(path, logging):
    if not os.path.exists(path):
        os.makedirs(path)
        logging.info(f"Created directory: {path}")
    else:
        logging.info(f"Directory already exists: {path}")


def read_csv_to_dict(file_path):
    data_dict = {}
    with open(file_path, "r") as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            variable = row.pop(
                ""
            )  # remove empty key from dictionary (first column in csv file)
            data_dict[variable] = {key: float(value) for key, value in row.items()}
    return data_dict


def load_reccap_mask():
    reccap_mask = xr.open_dataset(
        "./observations/RECCAP_AfricaSplit_MASK11_Mask_regridded.hadcm3bl_grid.nc"
    )
    regions = {
        1: "North_America",
        2: "South_America",
        3: "Europe",
        4: "Africa",  # combine North Africa (=4) and South Africa (+5)
        6: "North_Asia",
        7: "Central_Asia",
        8: "East_Asia",
        9: "South_Asia",
        10: "South_East_Asia",
        11: "Oceania",
    }
    return reccap_mask, regions


def _add_reccap_metrics(key, attributes, output_metrics):
    regions = {
        "North_America": "NAM",
        "South_America": "SAM",
        "Europe": "EUR",
        "Africa": "AFR",
        "North_Asia": "NAS",
        "Central_Asia": "CAS",
        "East_Asia": "EAS",
        "South_Asia": "SAS",
        "South_East_Asia": "SEA",
        "Oceania": "OCN",
    }

    cmip6_mean_values = read_csv_to_dict("./observations/stores_vs_fluxes_cmip6.csv")
    cmip6_error_values = read_csv_to_dict(
        "./observations/stores_vs_fluxes_cmip6_err.csv"
    )
    reccap_mean_values = read_csv_to_dict("./observations/stores_vs_fluxes_reccap.csv")
    reccap_error_values = read_csv_to_dict(
        "./observations/stores_vs_fluxes_reccap_err.csv"
    )

    for region, abbr in regions.items():
        if key == "RECCAP_VEG_C":
            var_name = f"RECCAP_{region}_sum_VEG_C"
            metric_realm = "global_carbon_stores"
            csv_var = "CVeg"
        elif key == "RECCAP_SOIL_C":
            var_name = f"RECCAP_{region}_sum_SOIL_C"
            metric_realm = "global_carbon_stores"
            csv_var = "CSoil"
        elif key == "RECCAP_GPP":
            var_name = f"RECCAP_{region}_sum_GPP"
            metric_realm = "global_productivity_fluxes"
            csv_var = "GPP"
        elif key == "RECCAP_TAU":
            var_name = f"RECCAP_{region}_sum_RH"
            metric_realm = "global_productivity_fluxes"
            csv_var = "Tau"

        target_min = min(
            cmip6_mean_values[csv_var][region] - cmip6_error_values[csv_var][region],
            reccap_mean_values[csv_var][region] - reccap_error_values[csv_var][region],
        )
        target_max = max(
            cmip6_mean_values[csv_var][region] + cmip6_error_values[csv_var][region],
            reccap_mean_values[csv_var][region] + reccap_error_values[csv_var][region],
        )

        # add metric to output dictionary
        output_metrics[f"{csv_var}_{abbr}"] = {
            "var_name": var_name,
            "metric_realm": metric_realm,
            "target_min": target_min,
            "target_max": target_max,
            "weight": attributes["weight"],
        }
    return output_metrics


def get_table_metrics(selected_metrics):

    table_metrics = {
        "GPP": {
            "var_name": "global_sum_GPP",
            "metric_realm": "global_productivity_fluxes",
            "target_min": 110,
            "target_max": 130,
        },
        "CVeg": {
            "var_name": "global_sum_VEG_C",
            "metric_realm": "global_carbon_stores",
            "target_min": 500,
            "target_max": 600,
        },
        "CSoil": {
            "var_name": "global_sum_SOIL_C",
            "metric_realm": "global_carbon_stores",
            "target_min": 1000,
            "target_max": 1500,
        },
        "Tr30SN": {
            "var_name": "subtropical_trees_30S_30N",
            "metric_realm": "global_veg_fractions",
            "target_min": 0.27 - 0.05,
            "target_max": 0.27 + 0.05,
        },
        "Tr30-90N": {
            "var_name": "NH_trees_30N_90N",
            "metric_realm": "global_veg_fractions",
            "target_min": 0.20 - 0.05,
            "target_max": 0.20 + 0.05,
        },
        "BareSoil": {
            "var_name": "global_mean_bare_soil",
            "metric_realm": "global_veg_fractions",
            "target_min": 0.26 - 0.05,
            "target_max": 0.26 + 0.05,
        },
    }

    output_metrics = {}
    for key, attributes in selected_metrics.items():
        if key in table_metrics:
            output_metrics[key] = table_metrics[key].copy()
            output_metrics[key]["weight"] = attributes["weight"]
        elif key in ["RECCAP_VEG_C", "RECCAP_SOIL_C", "RECCAP_GPP", "RECCAP_TAU"]:
            output_metrics = _add_reccap_metrics(key, attributes, output_metrics)

    return output_metrics
