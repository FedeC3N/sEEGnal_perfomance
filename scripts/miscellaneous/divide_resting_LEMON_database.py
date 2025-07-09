#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The recordings of LEMON_database have the EO and EC conditions together in the same recording (EO-EC-EO-EC...)
Divide to obtain EC and EO individually

Created on Thu 22/11/2024

@author: Fede
"""
import os
import glob

import mne
from eeglabio.utils import export_mne_epochs

# Resting state conditions
config = {}
config['rs_EO_code'] = 'Stimulus/S200'
config['rs_EC_code'] = 'Stimulus/S210'

# Filter the subjects of interest
config['database'] = 'LEMON_database'
config['path'] = {'sourcedata': os.path.join('databases','LEMON_database','sourcedata')}
files = glob.glob(os.path.join(config['path']['sourcedata'],'sub*vhdr'))
files = [os.path.basename(current_files) for current_files in files]


for current_file in files:

    print('Working on ' + current_file)

    try:

        #######
        # EO
        #######
        # Read the file
        #dummy_current_file = os.path.join(config['path']['sourcedata'],current_file)
        #raw = mne.io.read_raw(dummy_current_file,preload=True)

        # Get the annotations
        #annotations = raw.annotations

        # Read the annotations and go one by one
        #events, event_dict = mne.events_from_annotations(raw,regexp=config['rs_EO_code'])

        # Get the epochs
        #epochs = mne.Epochs(
        #    raw,
        #    events=events,
        #    tmin=0,tmax=2,
        #    reject_by_annotation=False,
        #    baseline=None,
        #    preload=True)

        # Export
        #filename = current_file[:10] + '_ses-1' + '_task-EO_eeg.set'
        #out_filepath = os.path.join(config['path']['sourcedata'],filename)
        #export_mne_epochs(epochs,out_filepath)

        #######
        # EC
        #######
        # Read the file
        dummy_current_file = os.path.join(config['path']['sourcedata'],current_file)
        raw = mne.io.read_raw(dummy_current_file,preload=True)

        # Get the annotations
        annotations = raw.annotations


        # Read the annotations and go one by one
        events,event_dict = mne.events_from_annotations(raw,regexp=config['rs_EC_code'])

        # Get the epochs
        epochs = mne.Epochs(
            raw,
            events=events,
            tmin=0,tmax=2,
            reject_by_annotation=False,
            baseline=None,
            preload=True)

        # Export
        filename = current_file[:10] + '_ses-1' + '_task-EC_eeg.set'
        out_filepath = os.path.join(config['path']['sourcedata'],filename)
        export_mne_epochs(epochs,out_filepath)


    except:
        print('   ERROR')
