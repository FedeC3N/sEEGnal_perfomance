#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
script to run each step of the ETL proccess across all subjects of a databse

Created on Thu 17/05/2024

@author: Fede
"""

# Imports
import os
import sys

import scripts.s0a_run_ETL_on_databases.private.init as init
import scripts.s0a_run_ETL_on_databases.private.create_subject as create_subject
from scripts.s0a_run_ETL_on_databases.private.etl_wrapper import standardize,badchannel_detection,artifact_detection,final_qa,export_clean


# Add the path
sys.path.append(os.path.join('..','TSD','aimind.etl'))

#### PARAMETERS
# Select the database
selected_database = 'AI_Mind_database'

# What step to run: standardize, badchannel, artifact, final_qa, export_clean
run = [1,1,1,1,1]

# Init Lurtis configuration
config = init.init_lurtis()
config['database'] = selected_database

# Init the database
subjects_id, sessions_id = init.init_database(config)
subjects_id = subjects_id[0:2]

# Select the desired type of file
#config.pattern = "[0-9]-[0-9]{3}-[0-9]-[A-Z]_(1-EO|2-EC|3-EO|4-EC)_.*"
config['pattern'] = "[0-9]-[0-9]{3}-[0-9]-[A-Z]_(.*-EO)_.*"

# List of subjects with errors
errors = []

# Go through each subject
for current_subject_id in subjects_id:

    # Go through each session
    for current_session_id in sessions_id:

        # Create the subjects following AI-Mind protocol
        bids_path = create_subject.create_bids_path(config,current_subject_id,current_session_id)

        # If the file does not exist, pass
        if not(str(bids_path.__class__) =="<class 'mne_bids.path.BIDSPath'>"):
            if len(bids_path) == 0:
                continue

        print('Work with ' + current_subject_id + ' - ' + current_session_id)

        try:

            # Run the selected processes
            if run[0]:
                print('   Standardize')
                standardize(config, bids_path)
                print('')

            if run[1]:
                print('   Badchannel detection')
                badchannel_detection(config,bids_path)
                print('')

            if run[2]:
                print('   Artifact Detection')
                artifact_detection(config,bids_path)
                print('')

            if run[3]:
                print('   Final QA')
                final_qa(config,bids_path)
                print('')

            if run[4]:
                print('   Export Clean')
                export_clean(config,bids_path)

        except:
            print('ERROR')
            errors.append(current_subject_id + ' - ' + current_session_id)


        print()
        print()
        print()


# Finally print the errors to check
if len(errors):

    from pprint import pprint
    print('LIST OF ERRORS')
    pprint(errors)


