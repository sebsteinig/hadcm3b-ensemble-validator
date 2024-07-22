import os
import subprocess
import logging

# Define the file paths and connection details
home_dir = os.path.expanduser("~")
id_list = os.path.join(home_dir, "hadcm3b-ensemble-generator", "logs", "xpzn_generated_ids_20240720.log")
local_data_dir = "/export/silurian/array-01/wb19586/ensembles"
remote_data_dir = "/user/home/wb19586/dump2hold"
remote_host = "bc4login.acrc.bris.ac.uk"
remote_user = "wb19586"
suffixes = ["pc", "pd", "pj"]

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('sync_files.log'),
        logging.StreamHandler()
    ]
)

# Function to read IDs from a log file
def read_ids_from_log(file_path):
    with open(file_path, 'r') as file:
        return [line.strip() for line in file.readlines()]

# Function to create a directory if it doesn't exist
def create_directory(path):
    if not os.path.exists(path):
        os.makedirs(path)
        logging.info(f"Created directory: {path}")
    else:
        logging.info(f"Directory already exists: {path}")

# Function to download files using rsync with wildcard
def download_files(remote_dir, local_dir, id, suffix):
    remote_path = os.path.join(remote_dir, f"*{id}*{suffix}*")
    print(remote_path)
    # rsync_command = [
    #     "rsync", "-avz", "--ignore-existing",
    #     f"{remote_user}@{remote_host}:{remote_path}",
    #     local_dir
    # ]
    # try:
    #     subprocess.run(rsync_command, check=True)
    #     logging.info(f"Successfully synced files from {remote_path} to {local_dir}")
    # except subprocess.CalledProcessError as e:
    #     logging.error(f"Failed to sync files from {remote_path} to {local_dir}. Error: {e}")

# Main logic
def main():
    # Read the IDs from the log file
    ids = read_ids_from_log(id_list)
    logging.info(f"Read {len(ids)} IDs from log file")
    
    for id in ids:
        id_local_dir = os.path.join(local_data_dir, id)
        create_directory(id_local_dir)
        
        for suffix in suffixes:
            suffix_local_dir = os.path.join(id_local_dir, suffix)
            create_directory(suffix_local_dir)
            
            remote_dir = os.path.join(remote_data_dir, id, "datam")
            download_files(remote_dir, suffix_local_dir, id, suffix)

if __name__ == "__main__":
    logging.info("Script started")
    main()
    logging.info("Script finished")