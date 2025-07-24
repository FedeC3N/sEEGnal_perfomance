#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Convert clean raw EEGlab data to clean epoch EEGlab data

Created on 09/12/2024

@author: Fede
"""

# Imports
import os
import glob

import mne

from eeglabio.utils import export_mne_epochs

#### PARAMETERS
config = {'path':{}}
config['path']['root'] = os.path.join('databases','LEMON_database','derivatives','lemon','clean')
subjects = glob.glob(os.path.join(config['path']['root'],'sub*'))
subjects = [current_subject.split(os.path.sep)[-1] for current_subject in subjects]

# Go through each subject
for current_subject in subjects:

        # Get the files to modify
        files = glob.glob(os.path.join(config['path']['root'],current_subject,'ses-1','eeg','sub*set'))

        for current_file in files:

                # Read the data
                try:

                    print('Working with file ' + current_file)

                    # Read the original file
                    raw = mne.io.read_raw_eeglab(current_file,preload=True)

                    # Define as 4 seconds epochs
                    epoch_definition = {"length":4,"overlap":0,"padding":2}
                    raw = mne.make_fixed_length_epochs(raw,
                                                       duration=epoch_definition['length'],
                                                       preload=True,
                                                       reject_by_annotation=False,
                                                       overlap=epoch_definition['overlap']
                                                       )

                    # Export epochs
                    out_filepath = current_file[0:115] + '_new' +  current_file[115:]
                    # Export
                    export_mne_epochs(raw,str(out_filepath))
                except:
                    continue

