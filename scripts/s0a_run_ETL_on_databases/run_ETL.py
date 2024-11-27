#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
script to run each step of the ETL proccess across all subjects of a databse

Created on Thu 17/05/2024

@author: Fede
"""

# Imports
from sEEGnal.io.init_databases import init_database
from sEEGnal.io.bids import create_bids_path
from sEEGnal.standardize.standardize import standardize
from sEEGnal.preprocess.badchannel_detection import badchannel_detection
from sEEGnal.preprocess.artifact_detection import artifact_detection
from sEEGnal.io.export_clean import export_clean

#### PARAMETERS
# Select the database
database = 'LEMON_database'

# What step to run: standardize, badchannel, artifact, export_clean
run = [1,1,0,0]

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
                print('   Standardize', end='. ')
                results = standardize(config,current_file,bids_path)
                print(' Result ' + results['result'])

            if run[1]:
                print('   Badchannel detection', end='. ')
                results = badchannel_detection(config,bids_path)
                print(' Result ' + results['result'])

            if run[2]:
                print('   Artifact Detection', end='. ')
                results = artifact_detection(config, bids_path)
                print(' Result ' + results['result'])

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


