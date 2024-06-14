#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to plot raw and clean EEGs

Created on Wed 22/05/2024

@author: Fede
"""

# Imports
import os
import sys
sys.path.append(os.path.join('..','TSD','aimind.etl'))

import scripts.miscellaneous.private.init as init
import scripts.miscellaneous.private.create_subject as create_subject
import scripts.miscellaneous.private.plots as plots



#### PARAMETERS
# Select the database
selected_database = 'AI_Mind_database'
tasks = ['1EO', '2EC', '3EO', '4EC']



# Init Lurtis configuration
config = init.init_lurtis()
config['database'] = selected_database

# Init the database
subjects_id, sessions_id = init.init_database(config)

# Go through each subject
for current_subject_id in subjects_id:

    # Go through each session
    for current_session_id in sessions_id:

        for current_task in tasks:

            # Create the subjects following AI-Mind protocol
            bids_path = create_subject.create_bids_path(config,current_subject_id,current_session_id,current_task)

            # If the file does not exist, pass
            if str(bids_path.__class__) == "<class 'list'>":
                continue

            print('Work with ' + bids_path.basename)

            # Plots
            plots.plot_raw(bids_path)

            # Add to the final list
            while True:
                user_input = input('Valid user? (y/n): ')

                if user_input.lower() == 'y':

                    with open(os.path.join('.','data','AI_Mind_database','valid_subjects.txt'),'a') as f:
                        to_write = bids_path.basename + '\r'
                        f.write(to_write)

                    break
                elif user_input.lower() == 'n':
                    break
                else:
                    print('Type y or n')
