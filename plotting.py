import xarray as xr
import matplotlib.pyplot as plt
import os
import pandas as pd
import cftime
import json

# read job IDs
def read_ids_from_log(file):
    with open(file, 'r') as file:
        return [line.strip() for line in file.readlines()]
    
def get_model_params(json_file):
    with open(json_file, 'r') as file:
        data = json.load(file)
    # create a dictionary with ensemble_id as the key
    data_dict = {entry['ensemble_id']: entry for entry in data}
    return data_dict

def convert_time(ds):
    # time_units = ds['time'].attrs['units']
    time_units = "days since 1850-12-01 00:00:00"
    calendar = ds['time'].attrs.get('calendar', 'standard')
    time_values = ds['time'].values

    if isinstance(time_values[0], cftime.Datetime360Day):
        # Handle cftime datetime conversion to standard datetime
        time_values = cftime.num2date(time_values, units=time_units, calendar=calendar)
        time_values = [pd.Timestamp(t.strftime('%Y-%m-%d')) for t in time_values]
    else:
        # Handle standard datetime conversion
        reference_time_str = time_units.split('since')[1].strip()
        reference_time = pd.Timestamp(reference_time_str)
        time_values = pd.to_datetime(time_values, origin=reference_time, unit='D')

    return time_values


def plot_timeseries(metric, data_dir, experiment, output_dir, logging):

    ids = read_ids_from_log(f"./id_lists/{experiment}.log")

    # collect all data_vars from all files to determine the subplot grid
    all_vars = set()
    for id in ids:
        metric_file = os.path.join(data_dir, id, "processed", metric, f"{id}_global_productivity_fluxes.combined.nc")
        if os.path.isfile(metric_file):
            ds = xr.open_dataset(metric_file,decode_times=False)
            all_vars.update(v for v in ds.data_vars if v != 'time')

    # Create a plot with subplots for each data_var
    num_vars = len(all_vars)
    if num_vars == 0:
        logging.warning(f"No data variables found for metric {metric}. Skipping.")
        return

    fig, axes = plt.subplots(num_vars, 1, figsize=(10, 6 * num_vars))
    if num_vars == 1:
        axes = [axes]

    all_vars = sorted(all_vars)
    for ax, var in zip(axes, all_vars):
        for id in ids:
            metric_file = os.path.join(data_dir, id, "processed", metric, f"{id}_global_productivity_fluxes.combined.nc")
            if not os.path.isfile(metric_file):
                logging.warning(f"Metric file {metric_file} does not exist. Skipping.")
                continue

            ds = xr.open_dataset(metric_file, decode_times=False)
            if var in ds.data_vars:
                time = convert_time(ds)
                ax.plot(time, ds[var].values, label=f"{id}")

        # ax.set_xlabel("Time")
        ax.set_ylabel(ds[var].attrs['units'])
        ax.set_title(f"{var} timeseries for all ensemble IDs")
        
        # Place a single legend to the right of all subplots
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles, labels, loc='center left', bbox_to_anchor=(1, 0.5), fontsize='small')

    plt.tight_layout()
    plt.subplots_adjust(right=0.85)  # Adjust the right side to make space for the legend
    output_file = os.path.join(output_dir, f"{experiment}_{metric}_timeseries.png")
    plt.savefig(output_file)
    plt.close()
    
    logging.info(f"Saved plot for {metric} to {output_file}")

def plot_parameter_scatter(metric, data_dir, experiment, output_dir, logging):

    params = read_ids_from_log(f"./id_lists/{experiment}_updated_parameters.json")