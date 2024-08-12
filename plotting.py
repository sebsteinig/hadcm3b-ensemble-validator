import xarray as xr
import matplotlib.pyplot as plt
import os
import pandas as pd
import cftime
import numpy as np
import statsmodels.api as sm


def _get_variables(metric, model_params, data_dir):
    # collect all data_vars from all files to determine the subplot grid
    all_vars = []
    for id, params in model_params.items():
        metric_file = os.path.join(
            data_dir,
            id,
            "processed",
            metric,
            f"{id}_{metric}.combined.nc",
        )
        if os.path.isfile(metric_file):
            ds = xr.open_dataset(metric_file, decode_times=False)
            all_vars.extend(
                v for v in ds.data_vars if v != "time" and v not in all_vars
            )
    return all_vars


def _load_data(metric_file, var, logging):
    if not os.path.isfile(metric_file):
        logging.warning(f"Metric file {metric_file} does not exist. Skipping.")
        return None

    ds = xr.open_dataset(metric_file, decode_times=False)

    if ds.t.shape == ():
        logging.warning(f"Metric file {metric_file} has an empty time axis. Skipping.")
        return None

    if var in ds.data_vars:
        return ds

    return None


def _convert_time(ds):
    time_units = "days since 1850-12-01 00:00:00"
    calendar = ds["time"].attrs.get("calendar", "standard")
    time_values = ds["time"].values

    if isinstance(time_values[0], cftime.Datetime360Day):
        # Handle cftime datetime conversion to standard datetime
        time_values = cftime.num2date(time_values, units=time_units, calendar=calendar)
        time_values = [pd.Timestamp(t.strftime("%Y-%m-%d")) for t in time_values]
    else:
        # Handle standard datetime conversion
        reference_time_str = time_units.split("since")[1].strip()
        reference_time = pd.Timestamp(reference_time_str)
        time_values = pd.to_datetime(time_values, origin=reference_time, unit="D")

    return time_values


def add_observation_lines(ax, var, metric):
    if var == "global_sum_GPP":
        ax.axhline(
            y=110, color="k", linestyle="--", linewidth=2.0, label="target (min)"
        )
        ax.axhline(
            y=120, color="k", linestyle="-", linewidth=3.0, label="target (mean)"
        )
        ax.axhline(
            y=130, color="k", linestyle="--", linewidth=2.0, label="target (max)"
        )
    elif var == "global_sum_NPP":
        ax.axhline(y=50, color="k", linestyle="--", linewidth=2.0, label="target (min)")
        ax.axhline(y=60, color="k", linestyle="-", linewidth=3.0, label="target (mean)")
        ax.axhline(y=70, color="k", linestyle="--", linewidth=2.0, label="target (max)")
    elif var == "global_sum_VEG_C":
        ax.axhline(y=500, color="k", linestyle="-", linewidth=3.0, label="target (min)")
        ax.axhline(y=600, color="k", linestyle="-", linewidth=3.0, label="target (max)")
    elif var == "global_sum_SOIL_C":
        ax.axhline(y=1000, color="k", linestyle="-", linewidth=3.0, label="target (min)")
        ax.axhline(y=1500, color="k", linestyle="-", linewidth=3.0, label="target (max)")
    
    # get precomputed observation refrerence values
    if metric == "global_veg_fractions":
        ds_obs = xr.open_dataset("./observations/igbp.veg_fraction_metrics.nc", decode_times=False)
        ax.axhline(y=ds_obs[var] - 0.1 * ds_obs[var], color="k", linestyle="--", linewidth=2.0, label="obs (min)")
        ax.axhline(y=ds_obs[var], color="k", linestyle="-", linewidth=3.0, label="obs (mean)")
        ax.axhline(y=ds_obs[var] + 0.1 * ds_obs[var], color="k", linestyle="--", linewidth=2.0, label="obs (max)")


# adapted from https://stackoverflow.com/a/59756979/3565452
def _simple_regplot(
    x, y, n_std=2, n_pts=100, ax=None, scatter_kws=None, line_kws=None, ci_kws=None
):
    """Draw a regression line with error interval."""
    ax = plt.gca() if ax is None else ax

    # calculate best-fit line and interval
    x_fit = sm.add_constant(x)
    fit_results = sm.OLS(y, x_fit).fit()

    eval_x = sm.add_constant(np.linspace(np.min(x), np.max(x), n_pts))
    pred = fit_results.get_prediction(eval_x)

    # draw the fit line and error interval
    ci_kws = {} if ci_kws is None else ci_kws
    ax.fill_between(
        eval_x[:, 1],
        pred.predicted_mean - n_std * pred.se_mean,
        pred.predicted_mean + n_std * pred.se_mean,
        alpha=0.3,
        color='lightcoral',
    )

    ax.plot(eval_x[:, 1], pred.predicted_mean - n_std * pred.se_mean, 
            linestyle='--', color='lightcoral', **ci_kws)
    ax.plot(eval_x[:, 1], pred.predicted_mean + n_std * pred.se_mean, 
            linestyle='--', color='lightcoral', **ci_kws)
    
    line_kws = {} if line_kws is None else line_kws
    h = ax.plot(eval_x[:, 1], pred.predicted_mean, **line_kws)

    # draw the scatterplot
    scatter_kws = {} if scatter_kws is None else scatter_kws
    ax.scatter(x, y, **scatter_kws)

    return fit_results


def plot_timeseries(metric, model_params, data_dir, experiment, output_dir, logging):

    # ceate a plot with subplots for each data_var
    all_vars = _get_variables(metric, model_params, data_dir)
    num_vars = len(all_vars)
    if num_vars == 0:
        logging.warning(f"No data variables found for metric {metric}. Skipping.")
        return

    fig, axes = plt.subplots(num_vars, 1, figsize=(10, 6 * num_vars))
    if num_vars == 1:
        axes = [axes]

    for ax, var in zip(axes, all_vars):
        t_min = []
        t_max = []
        for id, params in model_params.items():
            metric_file = os.path.join(
                data_dir,
                id,
                "processed",
                metric,
                f"{id}_{metric}.combined.nc",
            )
            data = _load_data(metric_file, var, logging)
            if data is not None:
                time = _convert_time(data)
                units = data[var].attrs["units"]
                mask = np.isfinite(data[var].values)
                ax.plot(time[mask], data[var].values[mask], label=f"{id}")
                t_min.append(time[0])
                t_max.append(time[-1])

        ax.set_xlabel("model year")
        ax.set_ylabel(units)
        ax.set_title(f"{var}", fontsize=20, fontweight="bold")
        ax.set_xlim(min(t_min), max(t_max))

        # Place a single legend to the right of all subplots
        handles, labels = ax.get_legend_handles_labels()
        if len(handles) < 50:
            ax.legend(
                handles,
                labels,
                loc="center left",
                bbox_to_anchor=(1, 0.5),
                fontsize="small",
            )

        add_observation_lines(ax, var, metric)

    plt.tight_layout()
    if len(all_vars) == 2:
        plt.subplots_adjust(top=0.85)
    elif len(all_vars) == 4:
        plt.subplots_adjust(top=0.95)
    elif len(all_vars) == 11:
        plt.subplots_adjust(top=0.97)
    plt.suptitle(
        f"HadCM3BL-C / {experiment} / {metric}",
        fontsize=20,
        fontweight="regular",
        y=0.99,
    )
    if len(handles) < 50:
        plt.subplots_adjust(right=0.85)  # make space for the legend

    output_file = os.path.join(output_dir, f"{experiment}_{metric}_timeseries.pdf")
    plt.savefig(output_file)
    plt.close()

    logging.info(f"Saved plot for {metric} to {output_file}")


def plot_parameter_scatter(
    metric,
    model_params,
    data_dir,
    experiment,
    output_dir,
    logging,
    clim_start_year,
    clim_end_year,
):
    # Get all variables to be plotted
    all_vars = _get_variables(metric, model_params, data_dir)

    num_vars = len(all_vars)
    if num_vars == 0:
        logging.warning(f"No data variables found for metric {metric}. Skipping.")
        return

    # Number of parameters to plot
    param_keys = list(
        model_params[next(iter(model_params))].keys()
    )
    # clean pramater list
    param_keys.remove("ensemble_id") 
    if "TUPP" in param_keys:
        param_keys.remove("TUPP")
    num_params = len(param_keys)

    # Create subplots with one row per variable and one column per parameter
    fig, axes = plt.subplots(
        num_vars, num_params, figsize=(5 * num_params, 5 * num_vars), sharey="row"
    )
    if num_vars == 1:
        axes = [axes]
    if num_params == 1:
        axes = [[ax] for ax in axes]

    # Plot data for each variable and parameter
    for i, var in enumerate(all_vars):
        logging.info(
            f"Plotting parameter scatter for metric {metric} and variable {var}."
        )
        for j, param_key in enumerate(param_keys):
            ax = axes[i][j]
            x_values = []
            y_values = []

            for id, params in model_params.items():
                metric_file = os.path.join(
                    data_dir,
                    id,
                    "processed",
                    metric,
                    f"{id}_{metric}.combined.nc",
                )
                data = _load_data(metric_file, var, logging)
                if data is not None:
                    data["t"] = _convert_time(data)
                    y_value = float(
                        data[var]
                        .sel(t=slice(str(clim_start_year), str(clim_end_year)))
                        .mean()
                        .values
                    )
                    if np.isnan(y_value):
                        continue
                    x_value = params[param_key][0]  # focus on BL parameter
                    x_values.append(x_value)
                    y_values.append(y_value)
                    units = data[var].attrs["units"]

            fit = _simple_regplot(
                x_values,
                y_values,
                n_std=2,
                n_pts=100,
                ax=ax,
                scatter_kws={"color": "dodgerblue", "edgecolor": "black", "s": 50},
                line_kws={"color": "red", "linestyle": "-", "linewidth": 2},
                ci_kws=None,
            )
            ax.text(
                0.05,
                0.95,
                f"N={fit.nobs:.0f}\nSlope={fit.params[1]:.2f}\nIntercept={fit.params[0]:.2f}\n$R^2$={fit.rsquared**2:.2f}",
                transform=ax.transAxes,
                fontsize=10,
                verticalalignment="top",
                bbox=dict(facecolor="white", alpha=0.8),
            )

            ax.set_xlabel(param_key)
            if j == 0:
                ax.set_ylabel(f"{var} ({units})", fontsize=16, fontweight="bold")
            if i == 0:  # Only set title for the first row
                ax.set_title(param_key, fontsize=20, fontweight="bold")

            add_observation_lines(ax, var, metric)

    plt.tight_layout()
    if len(all_vars) == 2:
        plt.subplots_adjust(top=0.9)
    elif len(all_vars) == 4:
        plt.subplots_adjust(top=0.95)
    elif len(all_vars) == 11:
        plt.subplots_adjust(top=0.97)
    plt.suptitle(
        f"HadCM3BL-C / {experiment} / {metric} / {clim_start_year}-{clim_end_year} / BL",
        fontsize=20,
        fontweight="regular",
        y=0.99,
    )

    output_file = os.path.join(output_dir, f"{experiment}_{metric}_param_scatter.pdf")

    plt.savefig(output_file)
    plt.close(fig)  # Close the figure to avoid memory leaks

    logging.info(f"Saved plot for {metric} to {output_file}")

