#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
script to prepare the databases to run ETL

Created on Thu 17/05/2024

@author: Fede
"""




# Imports
import os
import re
import sys
import glob
import json



def init_database(database):
    """

    Call the correct init function

    """

    # Read the config dictionary
    with open('./config/etl_configuration.json','r') as file:
        config = json.load(file)
    config['database'] = database

    # Build the function name
    func = f"init_{database}"
    to_init = getattr(sys.modules[__name__],func)

    # Call the function
    config, files, sub, ses, task = to_init(config)

    return  config, files, sub, ses, task



def init_AI_Mind_database(config):

    # Folders to find the subjects
    config['path'] = {}
    config['path']['data_root'] = os.path.join('databases','AI_Mind_database')
    config['path']['sourcedata'] = os.path.join(config['path']['data_root'],'sourcedata')

    # Add extra info to the dictionary
    config['component_estimation']["channels_to_include"] = 'eeg'
    config['component_estimation']["channels_to_exclude"] = ['EMG1', 'EMG2', 'CLAV', 'VEOGL']
    config['badchannel_detection']["channels_to_include"] = 'eeg'
    config['badchannel_detection']["channels_to_exclude"] = ['EMG1', 'EMG2', 'CLAV', 'VEOGL']
    config['artifact_detection']["channels_to_include"] = 'eeg'
    config['artifact_detection']["channels_to_exclude"] = ['EMG1', 'EMG2', 'CLAV', 'VEOGL']

    # Filter the subjects of interest
    files = glob.glob(os.path.join(config['path']['sourcedata'],'*.cnt'))
    files = [os.path.basename(file) for file in files]

    # Extract the info
    pattern = "([0-9]-[0-9]{3})-([0-9])-[A-Z]_(1-EO|2-EC|3-EO|4-EC)_.*"
    matches = [re.findall(pattern,file) for file in files]
    matches,files = zip(*[(match,file) for match,file in zip(matches,files) if match])
    sub = [match[0][0] for match in matches]
    ses = [match[0][1] for match in matches]
    task = [match[0][2] for match in matches]

    return config, files, sub, ses, task
