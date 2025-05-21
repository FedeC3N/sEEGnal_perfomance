#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Helper script to estimate different information of interest

Created on Thu 21/05/2025

@author: Fede
"""

# Imports
import numpy

import mne_bids

import sEEGnal.io.bids as bids
import sEEGnal.tools.mnetools as aimind_mne
from sEEGnal.io.init_databases import init_database


def estimate_deviation(database):

    # Init the database
    config, files, sub, ses, task = init_database(database)

    # Create the output matrix
    if database == 'LEMON_database':
        deviation_all = numpy.empty((61, len(files)))
    else:
        deviation_all = numpy.empty((126,len(files)))

    # Go through each subject
    for current_index in range(len(files)):

        # current info
        current_file = files[current_index]
        current_sub = sub[current_index]
        current_ses = ses[current_index]
        current_task = task[current_index]

        # Create the subjects following AI-Mind protocol
        bids_path = bids.create_bids_path(config, current_file, current_sub, current_ses, current_task)

        # Parameters for loading EEG  recordings
        freq_limits = [config['badchannel_detection']['impossible_amplitude']['low_freq'],
                       config['badchannel_detection']['impossible_amplitude']['high_freq']]
        crop_seconds = [config['badchannel_detection']['crop_seconds']]
        channels_to_include = config['badchannel_detection']["channels_to_include"]
        channels_to_exclude = config['badchannel_detection']["channels_to_exclude"]

        # Load the raw data
        raw = aimind_mne.prepare_raw(
            config,
            bids_path,
            preload=True,
            channels_to_include=channels_to_include,
            channels_to_exclude=channels_to_exclude,
            freq_limits=freq_limits,
            crop_seconds=crop_seconds,
            badchannels_to_metadata=False,
            exclude_badchannels=False,
            set_annotations=False)

        # Select the current channels
        raw.pick(channels_to_include)

        # De-mean the channels
        raw_data = raw.get_data().copy()
        mean_per_channel = raw_data.mean(axis=1)
        raw_data_demean = raw_data - mean_per_channel[:, numpy.newaxis]

        # Estimate the average standard deviation of each channel
        raw_data_demean_std = raw_data_demean.std(axis=1)
        deviation_all[:,current_index] = raw_data_demean_std

    return deviation_all



def badchannel_types(database):

    # Init the database
    config, files, sub, ses, task = init_database(database)

    # Channel types
    ch_types = ['n/a','bad_high_impedance', 'bad_impossible_amplitude_badchannels', 'bad_power_spectrum',
                'bad_gel_bridge', 'bad_high_variance']

    # Go through each subject
    for current_index in range(len(files)):

        # current info
        current_file = files[current_index]
        current_sub = sub[current_index]
        current_ses = ses[current_index]
        current_task = task[current_index]

        # Create the subjects following AI-Mind protocol
        bids_path = bids.create_bids_path(config, current_file, current_sub, current_ses, current_task)

        # Builds the path to the file.
        tsv_file = bids.build_derivative(bids_path, 'channels.tsv')

        # Reads the contents of the channel file.
        ch_info = mne_bids.tsv_handler._from_tsv(tsv_file)

        # Print each type
        print('Working with sub ' + current_sub + ' ses ' + current_ses + ' task ' + current_task)
        for current_type in ch_types:

            # Get the indexes of the channels
            indexes = [i for i in range(len(ch_info['status_description'])) if
                       ch_info['status_description'][i] == current_type]

            # Get the channels
            current_channels = [ch_info['name'][i] for i in indexes]

            # Print
            if current_type == 'n/a':
                print(f"Good channels ({len(current_channels)}): {current_channels}")
            else:
                print(f"{current_type} ({len(current_channels)}): {current_channels}")
        print()