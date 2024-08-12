import os
import xarray as xr
import numpy as np
import glob


def _calculate_grid_cell_areas(ds):
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


def _load_data(data_dir, id, stream, logging):
    input_files_pattern = os.path.join(data_dir, id, f"{stream}/*.nc")
    input_files = glob.glob(input_files_pattern)
    
    if not input_files:
        logging.info(f"No files found for id {id} in {input_files_pattern}. Skipping processing.")
        return None

    ds = xr.open_mfdataset(input_files, combine='by_coords', decode_times=False).squeeze()

    return ds


def _clean_time(ds):
    if 't_1' in ds.dims:
        ds = ds.rename({'t_1': 't'})
    elif 't_2' in ds.dims:
        ds = ds.rename({'t_2': 't'})
    if 'time' in ds.dims:
        ds = ds.drop_dims('time')
    return ds
    

def _load_data_productivity(data_dir, id, stream, logging):
    input_files_pattern = os.path.join(data_dir, id, f"{stream}/*nv+.nc")
    input_files = sorted(glob.glob(input_files_pattern))
    
    if not input_files:
        logging.info(f"No files found for id {id} in {input_files_pattern}. Skipping processing.")
        return None
    
    # get template files
    pi_template = xr.open_dataset('./template_files/HadCM3BL-C.pi.mm.nc', decode_times=False)
    pi_template['GBMVegCarb_srf'].values = np.full_like(pi_template['GBMVegCarb_srf'].values, np.nan)
    pi_template['soilCarbon_mm_srf'].values = np.full_like(pi_template['soilCarbon_mm_srf'].values, np.nan)

    # Load all files iteratively and extract variables
    datasets_veg = []
    datasets_soil = []
    for file in input_files:
        ds = xr.open_dataset(file, decode_times=False)
        if 'soilCarbon_mm_srf' in ds.variables:
            ds_soil = _clean_time(ds['soilCarbon_mm_srf'])
        else:
            ds_soil = _clean_time(pi_template['soilCarbon_mm_srf'])
        ds_soil['t'] = ds['t']
        datasets_soil.append(ds_soil)

        if 'GBMVegCarb_srf' in ds.variables:
            ds_veg = _clean_time(ds['GBMVegCarb_srf'])
        else:
            ds_veg = _clean_time(pi_template['GBMVegCarb_srf'])
        ds_veg['t'] = ds['t']
        datasets_veg.append(ds_veg)

    ds_soil_combined = xr.concat(datasets_soil, dim='t').squeeze() if datasets_soil else None
    ds_veg_combined = xr.concat(datasets_veg, dim='t').squeeze() if datasets_veg else None

    if ds_soil_combined is not None and ds_veg_combined is not None:
        ds_combined = xr.merge([ds_soil_combined, ds_veg_combined])

    return ds_combined


def global_productivity_fluxes(data_dir, id, output_dir, logging):

    ds = _load_data(data_dir, id, "pj", logging)
    if ds is None:
        return
    
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
        global_sum = (ds[var_name] * area).sum(dim=['latitude', 'longitude'])
        ds_global_sum[f'global_sum_{flux}'] = global_sum * scaling_factor
        ds_global_sum[f'global_sum_{flux}'].attrs['units'] = 'PgC/yr'

    # write the results to a new NetCDF file
    output_file = f"{output_dir}/global_productivity_fluxes/{id}_global_productivity_fluxes.combined.nc"
    ds_global_sum.to_netcdf(output_file)

    logging.info(f"Global productivity timeseries for id {id} saved to {output_file}.")


def global_carbon_stores(data_dir, id, output_dir, logging):

    ds = _load_data_productivity(data_dir, id, "pi", logging)
    if ds is None:
        return
    
    variables = {
        'SOIL_C': 'soilCarbon_mm_srf',
        'VEG_C': 'GBMVegCarb_srf',
    }

    fluxes = ['VEG_C', 'SOIL_C']


    pfts = {
        0: 'BL',
        1: 'NL',
        2: 'C3',
        3: 'C4',
        4: 'Shrub'
    }

    ds_global_sum = xr.Dataset()
    ds_global_sum['time'] = ds['t']

    # calculate the global sum for each variable
    area = _calculate_grid_cell_areas(ds)
    # convert from kg C m-2 s-1 to PgC
    scaling_factor = 1e-12

    for flux in fluxes:
        var_name = variables[flux]
        
        # calculate the global sum over all latitudes and longitudes
        global_sum = (ds[var_name] * area).sum(dim=['latitude', 'longitude'])
        global_sum = xr.where(global_sum == 0, np.nan, global_sum)
        ds_global_sum[f'global_sum_{flux}'] = global_sum * scaling_factor
        ds_global_sum[f'global_sum_{flux}'].attrs['units'] = 'PgC'

    # write the results to a new NetCDF file
    output_file = f"{output_dir}/global_carbon_stores/{id}_global_carbon_stores.combined.nc"
    ds_global_sum.to_netcdf(output_file)

    logging.info(f"Global carbon stores timeseries for id {id} saved to {output_file}.")


def global_veg_fractions(data_dir, id, output_dir, logging):

    ds = _load_data_productivity(data_dir, id, "pi", logging)
    if ds is None:
        return
    
    variables = {
        'SOIL_C': 'soilCarbon_mm_srf',
        'VEG_C': 'GBMVegCarb_srf',
    }

    fluxes = ['VEG_C', 'SOIL_C']


    pfts = {
        0: 'BL',
        1: 'NL',
        2: 'C3',
        3: 'C4',
        4: 'Shrub'
    }

    ds_global_sum = xr.Dataset()
    ds_global_sum['time'] = ds['t']

    # calculate the global sum for each variable
    area = _calculate_grid_cell_areas(ds)
    # convert from kg C m-2 s-1 to PgC
    scaling_factor = 1e-12

    for flux in fluxes:
        var_name = variables[flux]
        
        # calculate the global sum over all latitudes and longitudes
        global_sum = (ds[var_name] * area).sum(dim=['latitude', 'longitude'])
        global_sum = xr.where(global_sum == 0, np.nan, global_sum)
        ds_global_sum[f'global_sum_{flux}'] = global_sum * scaling_factor
        ds_global_sum[f'global_sum_{flux}'].attrs['units'] = 'PgC'

    # write the results to a new NetCDF file
    output_file = f"{output_dir}/global_carbon_stores/{id}_global_carbon_stores.combined.nc"
    ds_global_sum.to_netcdf(output_file)

    logging.info(f"Global carbon stores timeseries for id {id} saved to {output_file}.")