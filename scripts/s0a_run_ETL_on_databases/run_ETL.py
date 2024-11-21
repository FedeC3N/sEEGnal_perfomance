#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
script to run each step of the ETL proccess across all subjects of a databse

Created on Thu 17/05/2024

@author: Fede
"""

# Imports
from etl.io.init_databases import init_database
from etl.io.bids import create_bids_path
from etl.standardize.standardize import standardize
from etl.preprocess.badchannel_detection import badchannel_detection
from etl.preprocess.artifact_detection import artifact_detection
from etl.io.export_clean import export_clean

#### PARAMETERS
# Select the database
database = 'AI_Mind_database'

# What step to run: standardize, badchannel, artifact, export_clean
run = [1,1,1,1]

# Init the database
config, files, sub, ses, task = init_database(database)

# List of subjects with errors
errors = []

# Go through each subject
for current_index in range(len(files)):

        # current info
        current_file = files[current_index]
        current_sub = sub[current_index]
        current_ses = ses[current_index]
        current_task = task[current_index]

        # Create the subjects following AI-Mind protocol
        bids_path = create_bids_path(config,current_file,current_sub,current_ses,current_task)

        print('Working with sub ' + current_sub + ' ses ' + current_ses + ' task ' + current_task)

        try:

            # Run the selected processes
            if run[0]:
                print('   Standardize')
                results = standardize(config,current_file,bids_path)

            if run[1]:
                print('   Badchannel detection')
                results = badchannel_detection(config,bids_path)

            if run[2]:
                print('   Artifact Detection')
                results = artifact_detection(config,bids_path)

            if run[3]:
                print('   Export clean')
                export_clean(config,bids_path)


        except:
            print('ERROR')
            errors.append(current_file)


        print()
        print()
        print()


# Finally print the errors to check
if len(errors):

    from pprint import pprint
    print('LIST OF ERRORS')
    pprint(errors)


