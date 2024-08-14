import os
import xarray as xr
import numpy as np
import glob

from common import load_reccap_mask

def _calculate_grid_cell_areas(ds):
    R = 6371e3  # Earth radius in meters

    # calculate the differences in latitude and longitude in radians
    dlat = np.deg2rad(np.abs(np.diff(ds["latitude"])))
    dlon = np.deg2rad(np.abs(np.diff(ds["longitude"])))

    # assume uniform grid spacing
    dlat = dlat[0]
    dlon = dlon[0]

    area = R**2 * np.cos(np.deg2rad(ds["latitude"])) * dlon * dlat
    # broadcast the area to match the shape of the dataset variables
    area = area * xr.ones_like(ds["longitude"])

    return area


def _regional_mean(ds, lon_bounds, lat_bounds, area):
    lat_values = ds["latitude"].values
    if not np.all(np.diff(lat_values) > 0):
        # if latitude is not increasing, reverse the latitude bounds
        lat_bounds = lat_bounds[::-1]

    region = ds.sel(
        longitude=slice(lon_bounds[0], lon_bounds[1]),
        latitude=slice(lat_bounds[0], lat_bounds[1]),
    )
    region_mean = region.weighted(
        area.sel(
            longitude=slice(lon_bounds[0], lon_bounds[1]),
            latitude=slice(lat_bounds[0], lat_bounds[1]),
        )
    ).mean(dim=["latitude", "longitude"])
    return region_mean


def _load_data(data_dir, id, stream, logging):
    input_files_pattern = os.path.join(data_dir, id, f"{stream}/*.nc")
    input_files = glob.glob(input_files_pattern)

    if not input_files:
        logging.info(
            f"No files found for id {id} in {input_files_pattern}. Skipping processing."
        )
        return None

    ds = xr.open_mfdataset(
        input_files, combine="by_coords", decode_times=False
    ).squeeze()

    return ds


def _clean_time(ds):
    if "t_1" in ds.dims:
        ds = ds.rename({"t_1": "t"})
    elif "t_2" in ds.dims:
        ds = ds.rename({"t_2": "t"})
    if "time" in ds.dims:
        ds = ds.drop_dims("time")
    return ds


def _load_data_pi(data_dir, id, stream, logging, variables):
    input_files_pattern_1 = os.path.join(data_dir, id, f"{stream}/*nv+.nc")
    input_files_pattern_2 = os.path.join(data_dir, id, f"{stream}/*nov+.nc")
    input_files = sorted(glob.glob(input_files_pattern_1) + glob.glob(input_files_pattern_2))

    if not input_files:
        logging.info(
            f"No files found for id {id} in {input_files_pattern_1}. Skipping processing."
        )
        return None

    # Load template files and initialize placeholders for all variables
    pi_template = xr.open_dataset(
        "./template_files/HadCM3BL-C.pi.mm.nc", decode_times=False
    )
    placeholder_templates = {}
    for var in variables:
        if var in pi_template.variables:
            placeholder_template = pi_template[var].copy()
            placeholder_template.values = np.full_like(
                placeholder_template.values, np.nan
            )
            placeholder_templates[var] = placeholder_template

    # Initialize dictionary to hold datasets for each variable
    combined_datasets = {var: [] for var in variables}

    # Load all files iteratively and extract variables
    for file in input_files:
        print(file)
        ds = xr.open_dataset(file, decode_times=False)

        for var in variables:
            if var in ds.variables:
                ds_var = _clean_time(ds[var])
            else:
                ds_var = _clean_time(placeholder_templates[var])
            ds_var["t"] = ds["t"] - 300.0  # account for November files
            combined_datasets[var].append(ds_var)

    # Combine all datasets for each variable along the 't' dimension
    combined_data = []
    for var in variables:
        if combined_datasets[var]:
            ds_combined = xr.concat(combined_datasets[var], dim="t").squeeze()
            combined_data.append(ds_combined)

    # Merge all combined datasets into one
    if combined_data:
        ds_final_combined = xr.merge(combined_data)
        return ds_final_combined

    return None
    

def global_productivity_fluxes(data_dir, id, output_dir, logging):

    ds = _load_data(data_dir, id, "pj", logging)
    if ds is None:
        return

    variables = {"GPP": "GPP_ym_srf", "NPP": "NPP_ym_srf", "RH": "soilResp_ym_srf"}
    # compute NEP
    ds["NEP"] = ds[variables["NPP"]] - ds[variables["RH"]]

    # list of variables to process
    fluxes = ["GPP", "NPP", "RH", "NEP"]

    # get the RECCAP mask for regional means
    reccap_mask, regions = load_reccap_mask()

    # create a new xarray.Dataset to store the results
    ds_global_sum = xr.Dataset()
    ds_global_sum["time"] = ds["t"]

    # calculate the global sum for each variable
    area = _calculate_grid_cell_areas(ds)
    # convert from kg C m-2 s-1 to PgC/yr
    scaling_factor = 3600 * 24 * 360 * 1e-12
    for flux in fluxes:
        if flux in variables:
            var_name = variables[flux]
        else:
            var_name = flux

        # calculate the global sum over all latitudes and longitudes
        global_sum = (ds[var_name] * area).sum(dim=["latitude", "longitude"])
        ds_global_sum[f"global_sum_{flux}"] = global_sum * scaling_factor
        ds_global_sum[f"global_sum_{flux}"].attrs["units"] = "PgC/yr"

        # also calculate the sum for each RECCAP region
        for region_id, region_name in regions.items():
            # create a mask for the current region
            if region_name == "Africa": # combine North Africa and South Africa
                region_mask = (reccap_mask["Region_Map"] == 4) | (reccap_mask["Region_Map"] == 5)
            else:
                region_mask = reccap_mask["Region_Map"] == region_id
            masked_area = area.where(region_mask)
            regional_sum = (ds[var_name] * masked_area).sum(dim=["latitude", "longitude"])
            
            # Store the regional mean in the dataset
            ds_global_sum[f"RECCAP_{region_name}_sum_{flux}"] = regional_sum * scaling_factor
            ds_global_sum[f"RECCAP_{region_name}_sum_{flux}"].attrs["units"] = "PgC/yr"

    # write the results to a new NetCDF file
    output_file = f"{output_dir}/global_productivity_fluxes/{id}_global_productivity_fluxes.combined.nc"
    ds_global_sum.to_netcdf(output_file)

    logging.info(f"Global productivity timeseries for id {id} saved to {output_file}.")


def global_carbon_stores(data_dir, id, output_dir, logging):

    vars_to_load = ["soilCarbon_mm_srf", "GBMVegCarb_srf"]
    ds = _load_data_pi(data_dir, id, "pi", logging, vars_to_load)
    if ds is None:
        return

    variables = {
        "SOIL_C": "soilCarbon_mm_srf",
        "VEG_C": "GBMVegCarb_srf",
    }

    fluxes = ["VEG_C", "SOIL_C"]

    # get the RECCAP mask for regional means
    reccap_mask, regions = load_reccap_mask()

    ds_global_sum = xr.Dataset()
    ds_global_sum["time"] = ds["t"]

    # calculate the global sum for each variable
    area = _calculate_grid_cell_areas(ds)
    # convert from kg C m-2 s-1 to PgC
    scaling_factor = 1e-12

    for flux in fluxes:
        var_name = variables[flux]

        # calculate the global sum over all latitudes and longitudes
        global_sum = (ds[var_name] * area).sum(dim=["latitude", "longitude"])
        global_sum = xr.where(global_sum == 0, np.nan, global_sum)
        ds_global_sum[f"global_sum_{flux}"] = global_sum * scaling_factor
        ds_global_sum[f"global_sum_{flux}"].attrs["units"] = "PgC"

        # also calculate the sum for each RECCAP region
        for region_id, region_name in regions.items():
            # create a mask for the current region
            if region_name == "Africa": # combine North Africa and South Africa
                region_mask = (reccap_mask["Region_Map"] == 4) | (reccap_mask["Region_Map"] == 5)
            else:
                region_mask = reccap_mask["Region_Map"] == region_id
            masked_area = area.where(region_mask)
            regional_sum = (ds[var_name] * masked_area).sum(dim=["latitude", "longitude"])
            regional_sum = xr.where(regional_sum == 0, np.nan, regional_sum)

            # Store the regional mean in the dataset
            ds_global_sum[f"RECCAP_{region_name}_sum_{flux}"] = regional_sum * scaling_factor
            ds_global_sum[f"RECCAP_{region_name}_sum_{flux}"].attrs["units"] = "PgC"

    # write the results to a new NetCDF file
    output_file = (
        f"{output_dir}/global_carbon_stores/{id}_global_carbon_stores.combined.nc"
    )
    ds_global_sum.to_netcdf(output_file)

    logging.info(f"Global carbon stores timeseries for id {id} saved to {output_file}.")


def _calculate_veg_metrics(ds, var_name, surface_name):

    ds_global_mean = xr.Dataset()
    ds_global_mean["time"] = ds["t"]

    # calculate the global sum for each variable
    area = _calculate_grid_cell_areas(ds)

    # calculate the global mean for each PFT
    pfts = {0: "BL", 1: "NL", 2: "C3", 3: "C4", 4: "shrub", 7: "bare_soil"}
    for pft_id, pft_name in pfts.items():
        data_2d = ds[var_name].isel({surface_name: pft_id})
        global_mean = data_2d.weighted(area).mean(dim=["latitude", "longitude"])
        
        ds_global_mean[f"global_mean_{pft_name}"] = global_mean
        # also store the 2D data for each PFT
        ds_global_mean[f"{pft_name}_2D"] = data_2d

    # some aggregates
    ds_global_mean[f"global_mean_trees"] = (
        ds_global_mean[f"global_mean_BL"] + ds_global_mean[f"global_mean_NL"]
    )
    ds_global_mean[f"global_mean_grass"] = (
        ds_global_mean[f"global_mean_C3"] + ds_global_mean[f"global_mean_C4"]
    )

    # some regional means
    # amazon trees
    lon_bounds = [290, 320]  # 290°E to 320°E
    lat_bounds = [-15, 5]  # 15°S to 5°N
    amazon_mean = _regional_mean(ds, lon_bounds, lat_bounds, area)
    amazon_trees = (
        amazon_mean[var_name].isel({surface_name: [0, 1]}).sum(dim=surface_name)
    )
    ds_global_mean[f"amazon_trees"] = amazon_trees

    # subtropical trees
    lon_bounds = [0, 366]
    lat_bounds = [-30, 30]
    subtropics_mean = _regional_mean(ds, lon_bounds, lat_bounds, area)
    subtropics_trees = (
        subtropics_mean[var_name].isel({surface_name: [0, 1]}).sum(dim=surface_name)
    )
    ds_global_mean[f"subtropical_trees_30S_30N"] = subtropics_trees

    # NH trees
    lon_bounds = [0, 360]
    lat_bounds = [30, 60]
    nh_mean = _regional_mean(ds, lon_bounds, lat_bounds, area)
    nh_trees = nh_mean[var_name].isel({surface_name: [0, 1]}).sum(dim=surface_name)
    ds_global_mean[f"NH_trees_30N_90N"] = nh_trees

    for var in ds_global_mean.data_vars:
        ds_global_mean[var].attrs["units"] = "fraction"

    return ds_global_mean


def global_veg_fractions_model(data_dir, id, output_dir, logging):
    # load and process model data
    vars_to_load = ["fracPFTs_mm_srf"]
    ds = _load_data_pi(data_dir, id, "pi", logging, vars_to_load)
    if ds is None:
        return

    ds_global_mean = _calculate_veg_metrics(ds, vars_to_load[0], "pseudo_2")

    # write the results to a new NetCDF file
    output_file = (
        f"{output_dir}/global_veg_fractions/{id}_global_veg_fractions.combined.nc"
    )
    ds_global_mean.to_netcdf(output_file)
    logging.info(
        f"Global vegetation fraction timeseries for id {id} saved to {output_file}."
    )


def global_veg_fractions_obs(obs_file, logging):
    ds_obs = xr.open_dataset(obs_file, decode_times=False).squeeze()
    if ds_obs is None:
        return

    ds_obs_mean = _calculate_veg_metrics(ds_obs, "fracPFTs_snp_srf", "pseudo")

    # write the results to a new NetCDF file
    output_file_obs = f"./observations/igbp.veg_fraction_metrics.nc"
    ds_obs_mean.to_netcdf(output_file_obs)
    logging.info(
        f"Global vegetation reference metrics for IGBP saved to {output_file_obs}."
    )
