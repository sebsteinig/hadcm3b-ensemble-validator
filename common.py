import os
import json

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
