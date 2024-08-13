import os
import json
import xarray as xr

def parse_jobs(filename):
    """
    Parses a JSON file and returns a dictionary with ensemble_id as keys
    and their corresponding model parameters as values.
    
    Args:
        filename (str): The path to the JSON file.
        
    Returns:
        dict: A dictionary with ensemble_id as keys and parameter entries as values.
    """
    with open(filename, 'r') as file:
        data = json.load(file)
    
    result = {}
    for entry in data:
        ensemble_id = entry.get('ensemble_id')
        if ensemble_id:
            result[ensemble_id] = entry
    
    return result


def create_directory(path, logging):
    if not os.path.exists(path):
        os.makedirs(path)
        logging.info(f"Created directory: {path}")
    else:
        logging.info(f"Directory already exists: {path}")


def load_reccap_mask():
    reccap_mask = xr.open_dataset("./observations/RECCAP_AfricaSplit_MASK11_Mask_regridded.hadcm3bl_grid.nc")
    regions = {
        1: "North_America",
        2: "South_America",
        3: "Europe",
        4: "Africa", # combine North Africa (=4) and South Africa (+5)
        6: "North_Asia",
        7: "Central_Asia",
        8: "East_Asia",
        9: "South_Asia",
        10: "South_East_Asia",
        11: "Oceania",
    }
    return reccap_mask, regions
