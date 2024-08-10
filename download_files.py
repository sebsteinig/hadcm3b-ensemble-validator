import os
import subprocess
import logging
import argparse

# UM output streams to download
streams = ["pi", "pj", "da"]
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


# read job IDs
def read_ids_from_log(file_path):
    with open(file_path, 'r') as file:
        return [line.strip() for line in file.readlines()]


def create_directory(path):
    if not os.path.exists(path):
        os.makedirs(path)
        logging.info(f"Created directory: {path}")
    else:
        logging.info(f"Directory already exists: {path}")


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

        
def main(experiment):
    id_list = f"./id_lists/{experiment}.log"
    ids = read_ids_from_log(id_list)
    logging.info(f"Read {len(ids)} IDs from log file")
    
    for id in ids:
        id_local_dir = os.path.join(local_data_dir, id)
        create_directory(id_local_dir)
        
        for stream in streams:
            stream_local_dir = os.path.join(id_local_dir, stream)
            create_directory(stream_local_dir)
            
            remote_dir = os.path.join(remote_data_dir, id, "datam")
            #remote_dir = os.path.join(remote_data_dir, id)

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


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Download and process UM ensemble output files.')
    parser.add_argument('experiments', nargs='+', help='List of experiment names to process')
    
    args = parser.parse_args()
    
    for experiment in args.experiments:
        logging.info(f"Ensemble download script started for experiment: {experiment}")
        main(experiment)
        logging.info(f"Ensemble download script finished for experiment: {experiment}")
