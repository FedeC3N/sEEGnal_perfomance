#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Helper script to estimate different information of interest

Created on Thu 21/05/2025

@author: Fede
"""

# Imports
import numpy

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
        freq_limits = [config['badchannel_detection']['impossible_amplitude_detection_low_freq'],
                       config['badchannel_detection']['impossible_amplitude_detection_high_freq']]
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