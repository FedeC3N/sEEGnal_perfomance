#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
script to run each step of the ETL proccess across all subjects of a databse

Created on Thu 17/05/2024

@author: Fede
"""

# Imports
import etl.init.init as init
from etl.init.create_bids_path import create_bids_path
from etl.standardize.standardize import standardize
from etl.preprocess.badchannel_detection import badchannel_detection

#### PARAMETERS
# Select the database
database = 'AI_Mind_database'

# What step to run: standardize, badchannel, artifact, final_qa
run = [1,1,0,0]

# Init the database
config, files, sub, ses, task = init.init_database(database)

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
                print('')

            if run[1]:
                print('   Badchannel detection')
                results = badchannel_detection(config,bids_path)
                print('')

            if run[2]:
                print('   Artifact Detection')
                artifact_detection(config,bids_path)
                print('')

            if run[3]:
                print('   Final QA')
                final_qa(config,bids_path)
                print('')


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


