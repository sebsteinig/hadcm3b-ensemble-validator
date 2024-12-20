import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.lines as mlines
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import os
import pandas as pd
import cftime
import numpy as np
import statsmodels.api as sm
import csv
import cartopy.crs as ccrs

from common import load_reccap_mask, read_csv_to_dict
from tqdm import tqdm

highlight_colors = [
    "black", 
    "tab:red", 
    "tab:green", 
    "tab:orange", 
    "tab:purple", 
    "tab:brown", 
    "tab:pink", 
    "tab:gray", 
    "tab:olive", 
    "tab:cyan", 
    "tab:blue", 
    "yellow", 
    "lightblue", 
    "darkgreen", 
    "magenta", 
    "lightgray", 
    "lime", 
    "darkred", 
    "darkorange", 
    "teal", 
    "gold",
    "lightcoral",
    "darkviolet",
    "deepskyblue"
]

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
                v
                for v in ds.data_vars
                # if v != "time" and v not in all_vars and not v.startswith("RECCAP")
                if v != "time" and v not in all_vars
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
            y=102, color="k", linestyle="--", linewidth=2.0, label="target (min)"
        )
        ax.axhline(
            y=123, color="k", linestyle="-", linewidth=3.0, label="target (mean)"
        )
        ax.axhline(
            y=135, color="k", linestyle="--", linewidth=2.0, label="target (max)"
        )
    elif var == "global_sum_NPP":
        ax.axhline(y=50, color="k", linestyle="--", linewidth=2.0, label="target (min)")
        ax.axhline(y=60, color="k", linestyle="-", linewidth=3.0, label="target (mean)")
        ax.axhline(y=70, color="k", linestyle="--", linewidth=2.0, label="target (max)")
    elif var == "global_sum_VEG_C":
        ax.axhline(y=380, color="k", linestyle="-", linewidth=3.0, label="target (min)")
        ax.axhline(y=536, color="k", linestyle="-", linewidth=3.0, label="target (max)")
    elif var == "global_sum_SOIL_C":
        ax.axhline(
            y=1000, color="k", linestyle="-", linewidth=3.0, label="target (min)"
        )
        ax.axhline(
            y=1500, color="k", linestyle="-", linewidth=3.0, label="target (max)"
        )
    # get precomputed observation refrerence values (assuming +/- 5% error)
    if metric == "global_veg_fractions":
        ds_obs = xr.open_dataset(
            "./observations/igbp.veg_fraction_metrics.nc", decode_times=False
        )
        ax.axhline(
            y=ds_obs[var] - 0.05,
            color="k",
            linestyle="--",
            linewidth=2.0,
            label="obs (min)",
        )
        ax.axhline(
            y=ds_obs[var], color="k", linestyle="-", linewidth=3.0, label="obs (mean)"
        )
        ax.axhline(
            y=ds_obs[var] + 0.05,
            color="k",
            linestyle="--",
            linewidth=2.0,
            label="obs (max)",
        )


# adapted from https://stackoverflow.com/a/59756979/3565452
def _simple_regplot(
    x, y, n_std=2, n_pts=100, ax=None, scatter_kws=None, line_kws=None, ci_kws=None, highlight_ids=None, x_values_highlight=None, y_values_highlight=None, show_legend = True
):
    """Draw a regression line with error interval."""
    ax = plt.gca() if ax is None else ax

    # filter out nan values
    # x = np.array(x)
    # y = np.array(y)
    # mask = np.isfinite(x) & np.isfinite(y)
    # x = x[mask]
    # y = y[mask]

    # # calculate best-fit line and interval
    # x_fit = sm.add_constant(x)
    # fit_results = sm.OLS(y, x_fit).fit()

    # eval_x = sm.add_constant(np.linspace(np.min(x), np.max(x), n_pts))
    # pred = fit_results.get_prediction(eval_x)

    # # draw the fit line and error interval
    # ci_kws = {} if ci_kws is None else ci_kws
    # ax.fill_between(
    #     eval_x[:, 1],
    #     pred.predicted_mean - n_std * pred.se_mean,
    #     pred.predicted_mean + n_std * pred.se_mean,
    #     alpha=0.3,
    #     color="lightcoral",
    # )

    # ax.plot(
    #     eval_x[:, 1],
    #     pred.predicted_mean - n_std * pred.se_mean,
    #     linestyle="--",
    #     color="lightcoral",
    #     **ci_kws,
    # )
    # ax.plot(
    #     eval_x[:, 1],
    #     pred.predicted_mean + n_std * pred.se_mean,
    #     linestyle="--",
    #     color="lightcoral",
    #     **ci_kws,
    # )

    # line_kws = {} if line_kws is None else line_kws
    # h = ax.plot(eval_x[:, 1], pred.predicted_mean, **line_kws)

    # draw the scatterplot
    scatter_kws = {} if scatter_kws is None else scatter_kws
    ax.scatter(x, y, **scatter_kws)

    # highlight specific points
    if highlight_ids is not None:
        for j, (id, description) in enumerate(highlight_ids.items()):
            print(id)
            print(highlight_colors[j])
            ax.scatter(
                x_values_highlight[j],
                y_values_highlight[j],
                color=highlight_colors[j],
                edgecolor="black",  # black outline
                label=f"{id} ({description})",
                s=200,
                zorder = 100
            )
        # ax.legend(loc='upper right', fontsize='small')
        if show_legend :
            ax.legend(loc='upper left', fontsize='small', bbox_to_anchor=(1.05, 1), borderaxespad=0.)

    # return fit_results
    return None



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
        tile_y = 0.99
        plt.subplots_adjust(top=0.85)
    elif len(all_vars) == 4:
        plt.subplots_adjust(top=0.95)
        tile_y = 0.99
    elif len(all_vars) == 11:
        plt.subplots_adjust(top=0.97)
        tile_y = 0.99
    else:
        plt.subplots_adjust(top=0.995)
        tile_y = 0.999
    plt.suptitle(
        f"HadCM3BL-C / {experiment} / {metric}",
        fontsize=20,
        fontweight="regular",
        y=tile_y,
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
    highlight_ids=None,
):
    # Get all variables to be plotted
    all_vars = _get_variables(metric, model_params, data_dir)

    num_vars = len(all_vars)
    if num_vars == 0:
        logging.warning(f"No data variables found for metric {metric}. Skipping.")
        return

    # Number of parameters to plot
    param_keys = list(model_params[next(iter(model_params))].keys())
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
            x_values_highlight = []
            y_values_highlight = []

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
                    if x_value == -9999:
                        x_value = 0.5
                    x_values.append(x_value)
                    y_values.append(y_value)
                    units = data[var].attrs["units"]
                    if id in highlight_ids:
                        x_values_highlight.append(x_value)
                        y_values_highlight.append(y_value)

            fit = _simple_regplot(
                x_values,
                y_values,
                n_std=2,
                n_pts=100,
                ax=ax,
                scatter_kws={"color": "dodgerblue", "edgecolor": "black", "s": 50},
                line_kws={"color": "red", "linestyle": "-", "linewidth": 2},
                ci_kws=None,
                highlight_ids=highlight_ids,
                x_values_highlight=x_values_highlight,
                y_values_highlight=y_values_highlight,
            )
            if fit is not None:
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
        tile_y = 0.99
    elif len(all_vars) == 4:
        plt.subplots_adjust(top=0.95)
        tile_y = 0.99
    elif len(all_vars) == 11:
        plt.subplots_adjust(top=0.97)
        tile_y = 0.99
    else:
        plt.subplots_adjust(top=0.995)
        tile_y = 0.999
    plt.suptitle(
        f"HadCM3BL-C / {experiment} / {metric} / {clim_start_year}-{clim_end_year} / BL",
        fontsize=20,
        fontweight="regular",
        y=tile_y,
    )

    output_file = os.path.join(output_dir, f"{experiment}_{metric}_param_scatter.pdf")

    plt.savefig(output_file)
    plt.close(fig)  # Close the figure to avoid memory leaks

    logging.info(f"Saved plot for {metric} to {output_file}")


def _process_data(data, var, clim_start_year, clim_end_year):
    data["t"] = _convert_time(data)
    y_value = float(
        data[var].sel(t=slice(str(clim_start_year), str(clim_end_year))).mean().values
    )
    if np.isnan(y_value):
        return None
    return y_value


def _get_RECCAP_data(
    model_params, data_dir, region, realm, clim_start_year, clim_end_year, logging, highlight_ids=None
):
    x_values = []
    y_values = []
    x_values_highlight = []
    y_values_highlight = []

    for id, params in model_params.items():
        flux_file = os.path.join(
            data_dir,
            id,
            "processed",
            f"global_productivity_fluxes",
            f"{id}_global_productivity_fluxes.combined.nc",
        )
        store_file = os.path.join(
            data_dir,
            id,
            "processed",
            f"global_carbon_stores",
            f"{id}_global_carbon_stores.combined.nc",
        )
        flux_ds = _load_data(flux_file, f"{region}_sum_GPP", logging)
        store_ds = _load_data(store_file, f"{region}_sum_VEG_C", logging)
        if flux_ds is not None and store_ds is not None:
            if realm == "veg":
                flux_clim = _process_data(
                    flux_ds, f"{region}_sum_GPP", clim_start_year, clim_end_year
                )
                store_clim = _process_data(
                    store_ds, f"{region}_sum_VEG_C", clim_start_year, clim_end_year
                )
            elif realm == "soil":
                flux_clim = _process_data(
                    flux_ds, f"{region}_sum_RH", clim_start_year, clim_end_year
                )
                store_clim = _process_data(
                    store_ds, f"{region}_sum_SOIL_C", clim_start_year, clim_end_year
                )
                # tau defined as CS/soil_resp
                if flux_clim is not None and store_clim is not None:
                    flux_clim = store_clim / flux_clim

            if flux_clim is not None and store_clim is not None:
                x_values.append(flux_clim)
                y_values.append(store_clim)
            if highlight_ids is not None and id in highlight_ids:
                x_values_highlight.append(flux_clim)
                y_values_highlight.append(store_clim)
            else:
                continue

    return x_values, y_values, x_values_highlight[::-1], y_values_highlight[::-1]


def _add_rectangle(ax, x_mean, x_error, y_mean, y_error):
    # xalculate the rectangle boundaries
    x_min, x_max = x_mean - x_error, x_mean + x_error
    y_min, y_max = y_mean - y_error, y_mean + y_error

    # Add a gray rectangle
    rect = patches.Rectangle(
        (x_min, y_min),  # bottom left corner
        x_max - x_min,  # width
        y_max - y_min,  # height
        linewidth=1,
        edgecolor="none",
        facecolor="gray",
        alpha=0.5,
    )
    ax.add_patch(rect)


def _add_cross(ax, x_mean, x_error, y_mean, y_error):
    # add a black cross for the RECCAP range
    ax.add_line(
        mlines.Line2D(
            [x_mean - x_error, x_mean + x_error],
            [y_mean, y_mean],
            color="black",
            linewidth=3,
        )
    )
    ax.add_line(
        mlines.Line2D(
            [x_mean, x_mean],
            [y_mean - y_error, y_mean + y_error],
            color="black",
            linewidth=3,
        )
    )


def _add_legend(ax):
    # Retrieve existing legend handles and labels
    handles, labels = ax.get_legend_handles_labels()
    
    # Create custom legend handles
    rectangle_patch = patches.Patch(color="gray", alpha=0.5, label="CMIP6")
    cross_lines = mlines.Line2D([], [], color="black", linewidth=2, label="RECCAP2")
    
    # Add the custom handles to the existing ones
    handles.extend([rectangle_patch, cross_lines])
    labels.extend(["CMIP6", "RECCAP2"])
    
    # Update the legend with the new combined handles and labels
    ax.legend(handles=handles, labels=labels, fontsize="small")



def _draw_RECCAP_scatter(x_values, y_values, ax, region, realm, markersize, highlight_ids, x_values_highlight, y_values_highlight):
    # add ensemble data
    ax.scatter(
        x_values,
        y_values,
        color="dodgerblue",
        edgecolor="black",
        s=markersize,
        zorder=100,
    )

    # highlight specific points
    if highlight_ids is not None:
        for j, (id, description) in enumerate(highlight_ids.items()):
            print(f"{id} ({description})")
            print(x_values_highlight[j],y_values_highlight[j])
            ax.scatter(
                x_values_highlight[j],
                y_values_highlight[j],
                color=highlight_colors[j],
                edgecolor="black",  # black outline
                label=f"{id} ({description})",
                s=markersize*2,
                zorder = 100
            )
        ax.legend(loc='upper right', fontsize='small')

    # read in reference data from CSV
    cmip6_mean_values = read_csv_to_dict("./observations/stores_vs_fluxes_cmip6.csv")
    cmip6_error_values = read_csv_to_dict(
        "./observations/stores_vs_fluxes_cmip6_err.csv"
    )
    reccap_mean_values = read_csv_to_dict("./observations/stores_vs_fluxes_reccap.csv")
    reccap_error_values = read_csv_to_dict(
        "./observations/stores_vs_fluxes_reccap_err.csv"
    )

    if realm == "veg":
        # ax.set_xlabel("GPP (PgC yr$^{-1}$)")
        # ax.set_ylabel("Veg Carbon (PgC)")
        ax.set_xlabel("GPP (PgC yr$^{-1}$; Beer et al. 2010)")
        ax.set_ylabel("Veg Carbon (PgC); Erb et al. 2018")
        _add_rectangle(
            ax,
            cmip6_mean_values["GPP"][region],
            cmip6_error_values["GPP"][region],
            cmip6_mean_values["CVeg"][region],
            cmip6_error_values["CVeg"][region],
        )
        # _add_cross(
        #     ax,
        #     reccap_mean_values["GPP"][region],
        #     reccap_error_values["GPP"][region],
        #     reccap_mean_values["CVeg"][region],
        #     reccap_error_values["CVeg"][region],
        # )
        # instead of RECCAP totals, we will now use:
        # https://www.science.org/doi/10.1126/science.1184984
        # https://www.nature.com/articles/nature25138
        # as references (note much higher agreement with CMIP6 range)
        _add_cross(
            ax,
            123.0,
            8.0,
            458.0,
            78.0,
        )
    elif realm == "soil":
        ax.set_xlabel("tau (yrs)")
        ax.set_ylabel("Soil Carbon (PgC)")
        _add_rectangle(
            ax,
            cmip6_mean_values["Tau"][region],
            cmip6_error_values["Tau"][region],
            cmip6_mean_values["CSoil"][region],
            cmip6_error_values["CSoil"][region],
        )
        _add_cross(
            ax,
            reccap_mean_values["Tau"][region],
            reccap_error_values["Tau"][region],
            reccap_mean_values["CSoil"][region],
            reccap_error_values["CSoil"][region],
        )

    ax.set_title(region)

    _add_legend(ax)


def plot_RECCAP_stores_vs_fluxes(
    model_params,
    data_dir,
    experiment,
    output_dir,
    logging,
    clim_start_year,
    clim_end_year,
    highlight_ids
):

    # start with global plots of Veg carbon vs GPP and Soil carbon vs. tau (CS/soil_resp)
    fig, axes = plt.subplots(1, 2, figsize=(10, 5))

    print(highlight_ids)
    x_values, y_values, x_values_highlight, y_values_hightlight = _get_RECCAP_data(
        model_params, data_dir, "global", "veg", clim_start_year, clim_end_year, logging, highlight_ids
    )
    _draw_RECCAP_scatter(x_values, y_values, axes[0], "global", "veg", 50, highlight_ids, x_values_highlight, y_values_hightlight)

    x_values, y_values, x_values_highlight, y_values_hightlight = _get_RECCAP_data(
        model_params,
        data_dir,
        "global",
        "soil",
        clim_start_year,
        clim_end_year,
        logging,
        highlight_ids
    )
    _draw_RECCAP_scatter(x_values, y_values, axes[1], "global", "soil", 50, highlight_ids, x_values_highlight, y_values_hightlight)

    plt.tight_layout()
    plt.subplots_adjust(top=0.85)
    plt.suptitle(
        f"HadCM3BL-C / {experiment} / RECCAP_stores_vs_fluxes / {clim_start_year}-{clim_end_year} / global",
        fontsize=16,
        fontweight="regular",
        y=0.99,
    )

    output_file = os.path.join(output_dir, f"{experiment}_RECCAP_global_scatter.pdf")

    plt.savefig(output_file)
    plt.close(fig)

    logging.info(
        f"Saved RECCAP_stores_vs_fluxes global plot for {experiment} to {output_file}"
    )

    # now panel plots for each region, separated by realm
    realms = ["veg", "soil"]
    reccap_mask, regions = load_reccap_mask()

    for realm in realms:
        fig, axes = plt.subplots(4, 3, figsize=(15, 20))
        for i, region in enumerate(regions.values()):
            x_values, y_values, x_values_highlight, y_values_hightlight = _get_RECCAP_data(
                model_params,
                data_dir,
                f"RECCAP_{region}",
                realm,
                clim_start_year,
                clim_end_year,
                logging,
                highlight_ids
            )
            _draw_RECCAP_scatter(
                x_values, y_values, axes.flatten()[i], region, realm, 100, highlight_ids, x_values_highlight, y_values_hightlight
            )

        axes[3, 1].remove()
        axes[3, 2].remove()

        plt.tight_layout()
        plt.subplots_adjust(top=0.95)
        plt.suptitle(
            f"HadCM3BL-C / {experiment} / RECCAP_stores_vs_fluxes / {clim_start_year}-{clim_end_year} / regional",
            fontsize=16,
            fontweight="regular",
            y=0.99,
        )

        output_file = os.path.join(
            output_dir, f"{experiment}_RECCAP_regional_scatter_{realm}.pdf"
        )

        plt.savefig(output_file)
        plt.close(fig)

        logging.info(
            f"Saved RECCAP_stores_vs_fluxes regional plot for {experiment} and realm {realm} to {output_file}"
        )


def _calculate_skill_score(clim_value, target_min, target_max):
    target_mean = (target_min + target_max) / 2
    target_range = target_max - target_min
    # calculate the normalized distance from the target mean
    normalized_difference = abs(clim_value - target_mean) / target_range
    # calculate the skill score between 0 and 1
    # clim_value == target_mean: 1.0
    # clim_value == target_min or clim_value == target_max: 0.5
    # clim_value == target_mean - target_range or clim_value == target_max + target_range: 0.0
    # outside this range: 0.0
    skill_score = max(0.0, 1.0 - normalized_difference)
    return skill_score


def plot_overview_table(
    model_params,
    data_dir,
    experiment,
    output_dir,
    logging,
    table_metrics,
    clim_start_year,
    clim_end_year,
):
    rows = []
    hits = []

    # for id, params in model_params.items():
    for id, params in tqdm(model_params.items(), desc="Processing IDs"):
        # create a dictionary to store the data for each ensemble member
        row_data = {"ID": id}
        hit_data = {"ID": id}

        # get model parameters
        for key, value in params.items():
            # Get the first value for each key
            if isinstance(value, list):
                if key == "V_CRIT_ALPHA":
                    if value[0] == -9999:  
                        row_data["V_CRIT"] = np.nan
                    else:
                        row_data["V_CRIT"] = value[0]
                else:
                    row_data[key] = value[0]
            elif key == "ensemble_id":
                continue

        # check ensemble member performance against target values
        for metric_key, target in table_metrics.items():
            # load data
            metric_file = os.path.join(
                data_dir,
                id,
                "processed",
                target["metric_realm"],
                f"{id}_{target['metric_realm']}.combined.nc",
            )
            data = _load_data(metric_file, target["var_name"], logging)

            # if metric_key.startswith("rmse_BL"):
            #     print(target["var_name"])
            #     print(data)
            #     quit()
            # get climatology
            if data is not None:
                data["t"] = _convert_time(data)
                clim_value = float(
                    data[target["var_name"]]
                    .sel(t=slice(str(clim_start_year), str(clim_end_year)))
                    .mean()
                    .values
                )
            else:
                clim_value = None

            if metric_key.startswith("Tau_"):
                # load additional data for soil carbon to calculate tau = CS/soil_resp
                metric_file = os.path.join(
                    data_dir,
                    id,
                    "processed",
                    "global_carbon_stores",
                    f"{id}_global_carbon_stores.combined.nc",
                )
                cs_var_name = target["var_name"].replace("RH", "SOIL_C")
                data = _load_data(metric_file, cs_var_name, logging)
                if data is not None:
                    data["t"] = _convert_time(data)
                    cs_clim = float(
                        data[cs_var_name]
                        .sel(t=slice(str(clim_start_year), str(clim_end_year)))
                        .mean()
                        .values
                    )
                else:
                    cs_clim = None

                if clim_value is not None and cs_clim is not None:
                    clim_value = cs_clim / clim_value
                else:
                    clim_value = None

            if clim_value is not None:
                if metric_key.startswith("rmse_"):
                    # min RMSE of fractional data is 0 (highest skill), maximum is 1.0 (lowest skill)
                    skill_score = np.max([1.0 - 3.0 * clim_value, 0.0])
                    hit_data[metric_key] = skill_score
                    row_data[metric_key] = round(clim_value, 2)
                else:
                    # compare with target values
                    skill_score = _calculate_skill_score(
                        clim_value, target["target_min"], target["target_max"]
                    )
                    hit_data[metric_key] = skill_score
                    if clim_value <= 1.0:
                        row_data[metric_key] = round(clim_value, 2)
                    else:
                        row_data[metric_key] = round(clim_value, 1)
            else:
                row_data[metric_key] = np.nan
                hit_data[metric_key] = np.nan

        # calculate overall score based on hits/misses across all metrics
        overall_weight = sum([target["weight"] for target in table_metrics.values()])
        if all(np.isnan(row_data[metric_key]) for target in table_metrics.values()):
            row_data["overall_score"] = np.nan
            hit_data["overall_score"] = np.nan
        else:
            overall_score = (
                sum(
                    [
                        hit_data[target] * table_metrics[target]["weight"]
                        for target in table_metrics.keys()
                    ]
                )
                / overall_weight
            )

            row_data["overall_score"] = round(overall_score, 2)
            hit_data["overall_score"] = round(overall_score, 2)

        # append the ensemble member data as a new row to teh table
        rows.append(row_data)
        hits.append(hit_data)

    # Create a DataFrame from the list of rows
    df = pd.DataFrame(rows)
    df_hits = pd.DataFrame(hits)

    # sort the DataFrame by overall score
    df.sort_values(by="overall_score", ascending=False, inplace=True)

    # Set the ID column as the index (optional, you can skip this if you want to keep ID as a column)
    df.set_index("ID", inplace=False)
    df_hits.set_index("ID", inplace=True)

    # Save the DataFrame as a CSV file
    output_file_csv = os.path.join(output_dir, f"{experiment}_overview_table.csv")
    df.to_csv(output_file_csv, index=False)  # Save with ID as a column

    # Plot the table using ax.table for more customization
    fig, ax = plt.subplots(figsize=(10, 6))

    # Hide the axes
    ax.axis("tight")
    ax.axis("off")

    # Create a table at the current axes
    table = ax.table(
        cellText=df.values, colLabels=df.columns, cellLoc="center", loc="center"
    )

    # Customize the table
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1.2, 1.2)  # Adjust table size

    # Automatically adjust column widths based on content
    table.auto_set_column_width(col=list(range(len(df.columns))))

    print(df_hits)

    cmap = cm.RdYlGn
    norm = mcolors.Normalize(vmin=0, vmax=1.0)

    # Customize cell colors
    for (i, j), cell in table.get_celld().items():
        if j == 0:  # ID column
            cell.set_text_props(weight="bold")
            cell.set_facecolor("darkgray")
        elif j < 8:  # model parameter columns
            cell.set_facecolor("lightgray")
        if i == 0:  # header row
            cell.set_facecolor("dodgerblue")
            cell.set_text_props(color="white", weight="bold")
        # Check df_hits value for ID and column name of the current cell
        id_value = df.iloc[i - 1, 0]  # Adjust index by -1 for the row to match df
        column_name = df.columns[j]

        # color the cell based on the skill score
        if column_name in df_hits.columns and i > 0:
            skill_score = df_hits.loc[id_value, column_name]
            if np.isnan(skill_score) or np.isnan(df.iloc[i - 1, j]):
                cell.set_facecolor("white")
            else:
                cell_color = cmap(norm(skill_score))
                cell.set_facecolor(cell_color)
        elif column_name == "overall_score" and i > 0:
            skill_score = df_hits.loc[id_value, "overall_score"]
            if np.isnan(skill_score) or np.isnan(df.iloc[i - 1, j]):
                cell.set_facecolor("white")
            else:
                cell_color = cmap(norm(overall_score))
                cell.set_facecolor(cell_color)

    # save the table to disk
    output_file_png = os.path.join(output_dir, f"{experiment}_overview_table.pdf")
    plt.savefig(output_file_png, bbox_inches="tight", pad_inches=0.1)
    plt.close()

    logging.info(
        f"Saved overview table for {experiment} to {output_file_csv} and {output_file_png}"
    )


def plot_skill_score_scatter(
    model_params, experiment, output_dir, logging, clim_start_year, clim_end_year, highlight_ids=None,
):

    # Number of parameters to plot
    param_keys = list(model_params[next(iter(model_params))].keys())
    # clean pramater list
    param_keys.remove("ensemble_id")
    if "TUPP" in param_keys:
        param_keys.remove("TUPP")
    num_params = len(param_keys)

    # Create subplots with one row per variable and one column per parameter
    fig, axes = plt.subplots(1, num_params, figsize=(30, 5), sharey="row")

    for j, param_key in enumerate(param_keys):
        ax = axes[j]
        x_values = []
        y_values = []

        print(j)
        print(len(param_keys))
        if j == ( len(param_keys) - 1):
            show_legend = True
        else:
            show_legend = False
            
        # Read the CSV file
        skill_score_csv = os.path.join(output_dir, f"{experiment}_overview_table.csv")
        df = pd.read_csv(skill_score_csv)

        if param_key == "V_CRIT_ALPHA":
            search_key = "V_CRIT"
        else:
            search_key = param_key
        x_values = df[search_key].values
        y_values = df["overall_score"].values

        x_values_highlight = []
        y_values_highlight = []
        if highlight_ids is not None:
            for j, (id, description) in enumerate(highlight_ids.items()):
                id_index = df.index[df["ID"] == id].tolist()
                x_values_highlight.append(x_values[id_index[0]])
                y_values_highlight.append(y_values[id_index[0]])

        fit = _simple_regplot(
            x_values,
            y_values,
            n_std=2,
            n_pts=100,
            ax=ax,
            scatter_kws={"color": "dodgerblue", "edgecolor": "black", "s": 50},
            line_kws={"color": "red", "linestyle": "-", "linewidth": 2},
            ci_kws=None,
            highlight_ids=highlight_ids,
            x_values_highlight=x_values_highlight,
            y_values_highlight=y_values_highlight,
            show_legend = show_legend
        )
        if fit is not None:
            ax.text(
                0.05,
                0.05,
                f"N={fit.nobs:.0f}\nSlope={fit.params[1]:.2f}\nIntercept={fit.params[0]:.2f}\n$R^2$={fit.rsquared**2:.2f}",
                transform=ax.transAxes,
                fontsize=10,
                verticalalignment="bottom",
                bbox=dict(facecolor="white", alpha=0.8),
            )

        # ax.set_xlabel(param_key)
        if j == 0:
            ax.set_ylabel(f"overall skill score", fontsize=16, fontweight="bold")
        ax.set_title(param_key, fontsize=20, fontweight="bold")

    plt.subplots_adjust(top=0.8)
    plt.suptitle(
        f"HadCM3BL-C / {experiment} / ensemble skill score / {clim_start_year}-{clim_end_year} / BL",
        fontsize=20,
        fontweight="regular",
        y=0.98,
    )

    output_file = os.path.join(
        output_dir, f"{experiment}_overall_skill_core_param_scatter.pdf"
    )

    plt.savefig(output_file)
    plt.close(fig)

    logging.info(f"Saved scatter plot for overall skill score to {output_file}")


def find_geo_coords(ds):
    """
    Find the latitude and longitude variable names in a xarray dataset.

    Parameters
    ----------
    ds : xarray.Dataset
        Dataset to search for latitude and longitude variable names.

    Returns
    -------
    tuple
        A tuple containing the longitude and latitude variable names.

    Raises
    ------
    ValueError
        If the latitude or longitude variable names cannot be determined automatically.
    """
    possible_lat_names = ["latitude", "lat", "LAT", "Latitude"]
    possible_lon_names = ["longitude", "lon", "LON", "Longitude"]

    lat_name = None
    lon_name = None

    for var in ds.coords:
        if var in possible_lat_names:
            lat_name = var
        if var in possible_lon_names:
            lon_name = var

    if lat_name is None or lon_name is None:
        raise ValueError(
            "Could not automatically determine the latitude or longitude variable names."
        )

    return lon_name, lat_name


def add_cyclic_point(data, coord):
    """
    Add a cyclic (wrap-around) point to the data array and coordinate.
    """
    cyclic_data = xr.concat([data, data.isel({coord: 0})], dim=coord)
    cyclic_coord = np.append(data[coord], data[coord][0] + 360)
    cyclic_data[coord] = cyclic_coord
    return cyclic_data


def plot_filled_map(
    ax,
    data,
    type="contourf",
    cmap="viridis",
    levels=None,
    colorbar=True,
    extent=None,
    labels=True,
    title="",
    **kwargs,
):
    """
    Plot a filled map of a variable in an xarray dataset using standard matplotlib plotting.

    Parameters:
    - ax (matplotlib.axes.Axes): The axes to plot the map on.
    - data (xarray.DataArray): The data to plot.
    - type (str, optional): The type of plot to create. Default is "contourf".
    - cmap (matplotlib.colors.Colormap, optional): The colormap to use. Default is "viridis".
    - levels (array-like, optional): The explicit levels for contourf or pcolormesh plot. Default is None.
    - **kwargs: Additional keyword arguments to pass to the plot function.

    Returns:
    - p (matplotlib.contour.QuadContourSet or matplotlib.collections.QuadMesh): The plot object.

    """

    lon_name, lat_name = find_geo_coords(data)

    # Handle cyclic longitude if necessary
    if data[lon_name].max() > 180:
        data = add_cyclic_point(data, lon_name)

    if type == "contourf":
        p = ax.contourf(
            data[lon_name],
            data[lat_name],
            data,
            transform=ccrs.PlateCarree(),
            cmap=cmap,
            levels=levels,
            extend="both",
            **kwargs,
        )
    elif type == "pcolormesh":
        # Create a colormap based on levels for discrete color mapping
        if levels is not None:
            norm = mcolors.BoundaryNorm(levels, ncolors=cmap.N, clip=True)

        p = ax.pcolormesh(
            data[lon_name],
            data[lat_name],
            data,
            transform=ccrs.PlateCarree(),
            cmap=cmap,
            norm=norm,
            **kwargs,
        )

    # Set the geographic extent if specified
    if extent:
        ax.set_extent(extent, crs=ccrs.PlateCarree())

    # Add gridlines and labels
    gl = ax.gridlines(
        draw_labels=labels,
        linewidth=1,
        color="gray",
        alpha=0.5,
        linestyle="--",
        x_inline=False,
        y_inline=False,
    )
    gl.top_labels = False
    gl.xlabel_style = {"size": 8, "color": "black"}
    gl.ylabel_style = {"size": 8, "color": "black"}

    ax.set_title(title)

    if colorbar:
        cbar = plt.colorbar(p, ax=ax, orientation="vertical", pad=0.05, shrink=0.7)
        # cbar.set_label('variable name(K)')

    return p


def plot_PFT_maps(
    model_params,
    data_dir,
    experiment,
    output_dir,
    logging,
    clim_start_year,
    clim_end_year,
    highlight_ids=None,
):

    # ceate subplots, with 3 rows (trees, grass, bare soil) and one column for each ensemble member + one for the observations
    # pfts_to_plot = ["trees", "grass", "bare_soil"]
    pfts_to_plot = ["BL_2D", "NL_2D", "C3_2D", "C4_2D", "shrub_2D", "bare_soil_2D"]
    fig, axes = plt.subplots(
        len(pfts_to_plot),
        len(highlight_ids) + 1,
        figsize=(5 * (len(highlight_ids) + 1), 3.5 * len(pfts_to_plot)),
        subplot_kw={'projection': ccrs.Robinson()}
    )

    # plot the observations
    obs = xr.open_dataset('./observations/qrparm.veg.frac_igbp.pp.nc', decode_times=False)['fracPFTs_snp_srf'].squeeze()
    obs_remap = xr.open_dataset('./observations/qrparm.veg.frac_igbp.pp.hadcm3bl.nc', decode_times=False).squeeze()
    obs_metrics = xr.open_dataset("./observations/igbp.veg_fraction_metrics.nc", decode_times=False).squeeze()
    # pfts = {0: "BL", 1: "NL", 2: "C3", 3: "C4", 4: "shrub", 7: "bare_soil"}

    for i, pft in enumerate(pfts_to_plot):
        if pft == "trees":
            pft_frac = obs.isel(pseudo=0) + obs.isel(pseudo=1)
            obs_remap["trees"] = obs_remap['fracPFTs_snp_srf'].isel(pseudo=0) + obs_remap['fracPFTs_snp_srf'].isel(pseudo=1)
        elif pft == "grass":
            pft_frac = obs.isel(pseudo=2) + obs.isel(pseudo=3)
            obs_remap["grass"] = obs_remap['fracPFTs_snp_srf'].isel(pseudo=2) + obs_remap['fracPFTs_snp_srf'].isel(pseudo=3)
        elif pft == "bare_soil":
            pft_frac = obs.isel(pseudo=7)
            obs_remap["bare_soil"] = obs_remap['fracPFTs_snp_srf'].isel(pseudo=7)
        else:
            pft_frac = obs_metrics[pft]

        im = plot_filled_map(
            axes[i][0],
            pft_frac,
            type="pcolormesh",
            cmap=plt.cm.RdYlGn,
            levels=np.linspace(0, 1, 11),
            colorbar=False,
            extent=None,
            labels=False,
            title=f"{pft} / IGBP obs",
        )

    # plot the model data
    for j, (id, description) in enumerate(highlight_ids.items()):
        print(id)
        metric_file = os.path.join(
            data_dir,
            id,
            "processed",
            "global_veg_fractions",
            f"{id}_global_veg_fractions.combined.nc",
        )
        data = _load_data(metric_file, f"BL_2D", logging)
        if data is not None:
            data["t"] = _convert_time(data)
            clim = data.sel(t=slice(str(clim_start_year), str(clim_end_year))).mean("t")
            for i, pft in enumerate(pfts_to_plot):
                if pft == "trees":
                    pft_frac = clim["BL_2D"] + clim["NL_2D"]
                elif pft == "grass":
                    pft_frac = clim["C3_2D"] + clim["C4_2D"]
                elif pft == "bare_soil":
                    pft_frac = clim["bare_soil_2D"]
                else:
                    pft_frac = clim[pft]


                # rmse = np.sqrt(np.mean((pft_frac - obs_remap[pft])**2))

                im = plot_filled_map(
                    axes[i][j+1],
                    pft_frac,
                    type="pcolormesh",
                    cmap=plt.cm.RdYlGn,
                    levels=np.linspace(0, 1, 11),
                    colorbar=False,
                    extent=None,
                    labels=False,
                    # title=f"{pft} / {id} / {description} / rmse = {rmse:.2f}",
                    title=f"{pft} / {id} / {description}",
                )

    cbar_ax = fig.add_axes([0.2, 0.05, 0.6, 0.03])  # Adjust the axes dimensions as necessary
    cbar = fig.colorbar(im, cax=cbar_ax, orientation='horizontal', extend='neither', label='fractional coverage')

    plt.tight_layout()
    plt.subplots_adjust(top=0.93, bottom=0.1)
    plt.suptitle(
        f"HadCM3BL-C / {experiment} / PFT clims / {clim_start_year}-{clim_end_year}",
        fontsize=20,
        fontweight="regular",
        y=0.99,
    )

    # output_file = os.path.join(output_dir, f"{experiment}_PFT_maps.pdf")
    output_file = os.path.join(output_dir, f"{experiment}_PFT_maps.png")
    plt.savefig(output_file,dpi=200)
    plt.close()

    logging.info(f"Saved plot for PFT maps to {output_file}")
