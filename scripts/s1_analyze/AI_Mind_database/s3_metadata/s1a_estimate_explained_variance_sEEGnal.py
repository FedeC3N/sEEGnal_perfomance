#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Estimate the variance removed when cancelling some components

Created on Thu 20/05/2025

@author: Fede
"""

# Imports
import numpy

from sEEGnal.io import bids
from sEEGnal.tools import mnetools
from sEEGnal.io.init_databases import init_database



#### PARAMETERS
# Select the database
database = 'AI_Mind_database'

# Init the database
config, files, sub, ses, task = init_database(database)

# List of subjects with errors
errors = []

# Go through each subject
ratio_all = numpy.empty((len(files),1))
for current_index in range(len(files)):

        # current info
        current_file = files[current_index]
        current_sub = sub[current_index]
        current_ses = ses[current_index]
        current_task = task[current_index]

        # Create the subjects following AI-Mind protocol
        bids_path = bids.create_bids_path(config,current_file,current_sub,current_ses,current_task)

        print('Working with sub ' + current_sub + ' ses ' + current_ses + ' task ' + current_task)

        # Parameters for loading EEG recordings
        freq_limits = [config['component_estimation']['low_freq_limit'],
                       config['component_estimation']['high_freq_limit']]
        resample_frequency = config['artifact_detection']['resampled_frequency_estimate_component']
        epoch_definition = config['artifact_detection']['epoch_definition']
        channels_to_include = config['component_estimation']["channels_to_include"]
        channels_to_exclude = config['component_estimation']["channels_to_exclude"]

        # Load raw EEG
        raw = mnetools.prepare_raw(
            config,
            bids_path,
            preload=True,
            channels_to_include=channels_to_include,
            channels_to_exclude=channels_to_exclude,
            resample_frequency=resample_frequency,
            freq_limits=freq_limits,
            badchannels_to_metadata=True,
            exclude_badchannels=True,
            set_annotations=False)

        # Original variance
        raw_matrix = raw.get_data()
        original_variance = numpy.sum(numpy.var(raw_matrix, axis=1))

        # Get the mixing and unmixing matrix
        sobi = bids.read_sobi(bids_path,'sobi')
        mixing = sobi.mixing_matrix_.copy()
        unmixing = sobi.unmixing_matrix_.copy()

        # Keep brain and other components
        components_to_exclude = []
        if 'brain' in sobi.labels_.keys():
            components_to_exclude.append(sobi.labels_['eog'])
        if 'other' in sobi.labels_.keys():
            components_to_exclude.append(sobi.labels_['ecg'])
        components_to_exclude = sum(components_to_exclude, [])

        # If desired components, apply them.
        clean = raw.copy()
        if len(components_to_exclude) > 0:
            # Remove the eog components
            sobi.apply(clean,exclude=components_to_exclude)

        # Estimate the new variance
        clean_matrix = clean.get_data()
        clean_variance = numpy.sum(numpy.var(clean_matrix, axis=1))

        # Estimate the ratio
        ratio = clean_variance / original_variance
        ratio_all[current_index] = ratio

        print(f'  Original variance: {original_variance}')
        print(f'  Clean variance: {clean_variance}')
        print(f'  Ratio: {ratio}')
        components_to_exclude = [int(x) for x in components_to_exclude]
        print(f'  Components removed: {components_to_exclude}')
        print('  ')

print(ratio_all)