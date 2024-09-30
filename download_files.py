import os
import subprocess
import logging
import argparse
import glob
import re
import xarray as xr

from common import parse_jobs, create_directory

# UM output streams to download
streams = ["pi", "pj", "da"]
#streams = ["pi", "pj"]
# user settings
username = os.getlogin()
# storage on BRIDGE servers
local_data_dir = f"/export/silurian/array-01/{username}/ensembles"
# information on HPC used for model integration
remote_data_dir = f"/user/home/{username}/dump2hold"
remote_host = "bc4login.acrc.bris.ac.uk"
#remote_data_dir = f"/home/bridge/tw23150/umdata"
#remote_host = "silurian.ggy.bris.ac.uk"
remote_user = username

# set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('download_files_from_BC4.log'),
        logging.StreamHandler()
    ]
)


# get list of available HPC files
def list_remote_files(remote_dir, id, stream):
    if stream == "da":
        remote_path = os.path.join(remote_dir, f"{id}*#{stream}*c1+")
    else:
        remote_path = os.path.join(remote_dir, f"{id}*{stream}*")
    ssh_command = f"ssh {remote_user}@{remote_host} 'ls {remote_path}'"
    try:
        result = subprocess.check_output(ssh_command, shell=True, text=True)
        return result.split()
    except subprocess.CalledProcessError as e:
        logging.error(f"Failed to list files on remote server. Error: {e}")
        return []


# get list of local (i.e. previously processed) files
def list_converted_files(directory):
    return {os.path.splitext(filename)[0] for filename in os.listdir(directory) if filename.endswith('.nc')}


# download individual file from HPC
def download_file(remote_dir, local_dir, file):
    remote_path = os.path.join(remote_dir, file)
    rsync_command = [
        "rsync", "-av", "--ignore-existing",
        f"{remote_user}@{remote_host}:{remote_path}",
        local_dir
    ]
    try:
        subprocess.run(rsync_command, check=True)
        logging.info(f"Successfully synced file from {remote_path} to {local_dir}")
    except subprocess.CalledProcessError as e:
        logging.error(f"Failed to sync file from {remote_path} to {local_dir}. Error: {e}")


# convert individual file from PP to netCDF format
def convert_to_nc(file_path):
    nc_file_path = file_path + ".nc"
    try:
        conversion_command = ["um2nc", file_path]
        subprocess.run(conversion_command, check=True)
        delete_command = ["rm", "-f", file_path]
        subprocess.run(delete_command, check=True)
        logging.info(f"Successfully converted {file_path} to {nc_file_path}")
    except subprocess.CalledProcessError as e:
        logging.error(f"Failed to convert {file_path} to {nc_file_path}. Error: {e}")
        subprocess.run(["rm", "-f", file_path], check=True)
        subprocess.run(["rm", "-f", nc_file_path], check=True)


def _find_min_max_year(directory, pattern="*#pi*.nc"):
    files = glob.glob(os.path.join(directory, pattern))
    # regular expression to extract the year from the filename
    year_pattern = re.compile(r"#pi00000(\d{4})")

    years = []
    for file in files:
        filename = os.path.basename(file)
        match = year_pattern.search(filename)
        if match:
            year = int(match.group(1))
            years.append(year)

    if not years:
        return None, None

    min_year = min(years)
    max_year = max(years)

    return min_year, max_year


def _create_yearly_data(id, input_dir, output_dir, min_year, max_year, stream):
    create_directory(output_dir, logging)
    for year in range(min_year, max_year + 1):
        expected_filepath = f"{output_dir}/{id}a#{stream}00000{year}.nc"
        if not os.path.exists(expected_filepath):
            monthly_files_pattern = os.path.join(input_dir, stream, f"*{stream}00000{year}*.nc")
            monthly_files = glob.glob(monthly_files_pattern)
            if len(monthly_files) == 12:
                # Initialize an empty list to store the loaded datasets
                datasets = []
                
                # Loop through each monthly file and process it with xarray
                for f in monthly_files:
                    try:
                        # Load the dataset

                        variables_to_load = ['soilCarbon_mm_srf','fracPFTs_mm_srf', 'longitude', 'latitude', 'pseudo_2']
                        drop_variables = [var for var in xr.open_dataset(f, decode_times=False).variables if var not in variables_to_load]

                        ds = xr.open_mfdataset(f, combine='by_coords', decode_times=False, drop_variables=drop_variables).squeeze()
                        # Select the specific variable and reduce dimensions if needed
                        # ds_selected = ds.isel(surface=0, drop=True)
                        
                        # Append the selected data to the list of datasets
                        datasets.append(ds)
                        
                    except Exception as e:
                        logging.error(f"Failed to process {f}. Error: {e}")
                        continue
                
                if datasets:
                    # Combine all the datasets along the time dimension
                    average_ds = xr.concat(datasets, dim='time').mean(dim='time')
                    
                    try:
                        # Save the combined dataset to a new NetCDF file
                        average_ds.to_netcdf(expected_filepath)
                        logging.info(f"Successfully created {expected_filepath}.")
                    except Exception as e:
                        logging.error(f"Failed to save {expected_filepath}. Error: {e}")
                else:
                    logging.error(f"No datasets found for year {year}. Skipping file creation.") 


def main(experiment):
    model_params = parse_jobs(f"./id_lists/{experiment}_parameters.json")
    logging.info(f"Read {len(model_params)} IDs from log file")
    
    for id, params in model_params.items():
    # for id, params in reversed(list(model_params.items())):
        # if id != "xqaba":
        #     break
        id_local_dir = os.path.join(local_data_dir, id)
        create_directory(id_local_dir, logging)
        
        for stream in streams:
            stream_local_dir = os.path.join(id_local_dir, stream)
            create_directory(stream_local_dir, logging)
            
            remote_dir = os.path.join(remote_data_dir, id, "datam")

            # get list of available files from remote server
            remote_files = list_remote_files(remote_dir, id, stream)
            logging.info(f"Found {len(remote_files)} files on remote server for ID {id} and stream {stream}")

            # get list of already converted files
            converted_files = list_converted_files(stream_local_dir)
            
            # filter out files that have already been converted
            new_files = [file for file in remote_files if os.path.basename(file) not in converted_files]
            logging.info(f"Found {len(new_files)} new files to download and convert for ID {id} and stream {stream}")
            
            # download and convert only the new files
            if new_files:
                for file in new_files:
                    download_file(remote_dir, stream_local_dir, file)
                    convert_to_nc(os.path.join(stream_local_dir, os.path.basename(file)))

            # calculate annual means from the downloaded monthly files
            # if stream == "pi":
            #     min_year, max_year = _find_min_max_year(stream_local_dir, pattern="*#pi*.nc")
            #     if min_year and max_year:
            #         logging.info(f"Found min year {min_year} and max year {max_year} for ID {id} and stream {stream}")
            #         id_yearly_dir = os.path.join(local_data_dir, id, f"{stream}_ym")
            #         logging.info(f"Checking for yearly mean {stream} data in {id_yearly_dir}")
            #         _create_yearly_data(id, id_local_dir, id_yearly_dir, min_year, max_year, stream)                    
            #     else:
            #         logging.error(f"Failed to find min and max year for ID {id} and stream {stream}")



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Download and process UM ensemble output files.')
    parser.add_argument('experiments', nargs='+', help='List of experiment names to process')
    
    args = parser.parse_args()
    
    for experiment in args.experiments:
        logging.info(f"Ensemble download script started for experiment: {experiment}")
        main(experiment)
        logging.info(f"Ensemble download script finished for experiment: {experiment}")
