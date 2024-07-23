import os
import xarray as xr
import numpy as np
import glob


def calculate_grid_cell_areas(ds):
    R = 6371e3  # Earth radius in meters
    
    # calculate the differences in latitude and longitude in radians
    dlat = np.deg2rad(np.abs(np.diff(ds['latitude'])))
    dlon = np.deg2rad(np.abs(np.diff(ds['longitude'])))
    
    # assume uniform grid spacing
    dlat = dlat[0]
    dlon = dlon[0]
    
    area = R**2 * np.cos(np.deg2rad(ds['latitude'])) * dlon * dlat
    # broadcast the area to match the shape of the dataset variables
    area = area * xr.ones_like(ds['longitude'])
    
    return area


def global_productivity_fluxes(data_dir, id, output_dir, logging):
    input_files_pattern = os.path.join(data_dir, id, "pj/*.nc")
    input_files = glob.glob(input_files_pattern)
    
    if not input_files:
        logging.info(f"No files found for id {id} in {input_files_pattern}. Skipping processing.")
        return

    ds = xr.open_mfdataset(input_files, combine='by_coords', decode_times=False).squeeze()

    variables = {
        'GPP': 'GPP_ym_srf',
        'NPP': 'NPP_ym_srf',
        'RH': 'soilResp_ym_srf'
    }
    # compute NEP 
    ds['NEP'] = ds[variables['NPP']] - ds[variables['RH']]

    # list of variables to process
    fluxes = ['GPP', 'NPP', 'RH', 'NEP']

    # Create a new xarray.Dataset to store the results
    ds_global_sum = xr.Dataset()
    ds_global_sum['time'] = ds['t']

    # xalculate the global sum for each variable
    area = calculate_grid_cell_areas(ds)
    # convert from kg C m-2 s-1 to PgC/yr
    scaling_factor = 3600 * 24 * 360 * 1e-12
    for flux in fluxes:
        if flux in variables:
            var_name = variables[flux]
        else:
            var_name = flux
        
        # calculate the global sum over all latitudes and longitudes
        global_sum = (ds[var_name] * area).sum(dim=['latitude', 'longitude'])
        ds_global_sum[f'global_sum_{flux}'] = global_sum * scaling_factor

    # write the results to a new NetCDF file
    output_file = f"{output_dir}/{id}_global_productivity_fluxes.combined.nc"
    ds_global_sum.to_netcdf(output_file)

    logging.info(f"Global productivity timeseries for id {id} saved to {output_file}.")