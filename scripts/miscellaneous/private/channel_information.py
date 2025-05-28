#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Helper script to estimate different information of interest

Created on Thu 21/05/2025

@author: Fede
"""

# Imports
import numpy
import matplotlib.pyplot as plt

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
                'bad_gel_bridge', 'bad_high_deviation']

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



def plot_correlation_values(database):
    # Init the database
    config, files, sub, ses, task = init_database(database)

    # Create the output matrix
    if database == 'LEMON_database':
        deviation_all = numpy.empty((61, len(files)))
    else:
        deviation_all = numpy.empty((126, len(files)))

    # Go through each subject
    for current_index in range(len(files)):
        # current info
        current_file = files[current_index]
        current_sub = sub[current_index]
        current_ses = ses[current_index]
        current_task = task[current_index]

        # Create the subjects following AI-Mind protocol
        bids_path = bids.create_bids_path(config, current_file, current_sub, current_ses, current_task)

        # Read the ICA information
        sobi = bids.read_sobi(bids_path, 'sobi-badchannels')

        # Parameters for loading EEG  recordings
        freq_limits = [config['badchannel_detection']['gel_bridge']['low_freq'],
                       config['badchannel_detection']['gel_bridge']['high_freq']]
        channels_to_include = config['badchannel_detection']["channels_to_include"]
        channels_to_exclude = config['badchannel_detection']["channels_to_exclude"]
        crop_seconds = [config['badchannel_detection']["crop_seconds"]]
        epoch_definition = config['badchannel_detection']['gel_bridge']['epoch_definition']

        # Load the raw EEG
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
            set_annotations=False,
            epoch=epoch_definition)

        # Select the current channels
        raw.pick(channels_to_include)

        # If there is EOG or EKG, remove those components
        components_to_exclude = []
        if 'eog' in sobi.labels_.keys():
            components_to_exclude.append(sobi.labels_['eog'])
        if 'ecg' in sobi.labels_.keys():
            components_to_exclude.append(sobi.labels_['ecg'])
        components_to_exclude = sum(components_to_exclude, [])

        # If desired components, apply them.
        if len(components_to_exclude) > 0:
            # Remove the eog components
            sobi.apply(raw, exclude=components_to_exclude)

        # Get the data
        raw_data = raw.get_data()

        # Estimate the correlation matrix for each epoch
        correlation_coefficients = numpy.empty((raw_data.shape[1], raw_data.shape[1], raw_data.shape[0]))
        for i in range(raw_data.shape[0]):
            current_correlation = numpy.corrcoef(raw_data[i, :, :])
            mask = numpy.triu(numpy.ones_like(current_correlation, dtype=bool),k=1)
            current_correlation = numpy.where(mask,current_correlation,numpy.nan)
            correlation_coefficients[:, :, i] = current_correlation

        correlation_coefficients = correlation_coefficients.mean(axis=2)
        dummy = correlation_coefficients.reshape((1,-1))
        dummy = dummy[~numpy.isnan(dummy)]
        x = (current_index + 1 )*numpy.ones((len(dummy),1))
        plt.plot(x,dummy,'o')

    plt.show(block=True)




def compare_correlation_values(database):

        # Init the database
        config, files, sub, ses, task = init_database(database)

        print("SUBJECT - EPOCH - WHOLE LENGTH")

        # Go through each subject
        for current_index in range(len(files)):

            # current info
            current_file = files[current_index]
            current_sub = sub[current_index]
            current_ses = ses[current_index]
            current_task = task[current_index]

            # Create the subjects following AI-Mind protocol
            bids_path = bids.create_bids_path(config, current_file, current_sub, current_ses, current_task)

            # Create the output matrix
            badchannels_whole_length = []
            badchannels_epoched = []

            ###### EPOCH DATA

            # Read the ICA information
            sobi = bids.read_sobi(bids_path, 'sobi-badchannels')

            # Parameters for loading EEG  recordings
            freq_limits = [config['badchannel_detection']['gel_bridge']['low_freq'],
                           config['badchannel_detection']['gel_bridge']['high_freq']]
            channels_to_include = config['badchannel_detection']["channels_to_include"]
            channels_to_exclude = config['badchannel_detection']["channels_to_exclude"]
            crop_seconds = [config['badchannel_detection']["crop_seconds"]]
            epoch_definition = config['badchannel_detection']['gel_bridge']['epoch_definition']

            # Load the raw EEG
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
                set_annotations=False,
                epoch=epoch_definition)

            # If there is EOG or EKG, remove those components
            components_to_exclude = []
            if 'eog' in sobi.labels_.keys():
                components_to_exclude.append(sobi.labels_['eog'])
            if 'ecg' in sobi.labels_.keys():
                components_to_exclude.append(sobi.labels_['ecg'])
            components_to_exclude = sum(components_to_exclude, [])

            # If desired components, apply them.
            if len(components_to_exclude) > 0:
                # Remove the eog components
                sobi.apply(raw, exclude=components_to_exclude)

            # Get the data
            raw_data = raw.get_data().copy()

            # Estimate the correlation matrix for each epoch
            correlation_coefficients_mask = numpy.empty((raw_data.shape[1], raw_data.shape[1], raw_data.shape[0]))
            for i in range(raw_data.shape[0]):
                # I'm going to save the values above the threshold
                current_correlation = numpy.corrcoef(raw_data[i, :, :])
                current_correlation = numpy.triu(current_correlation, k=1)
                correlation_coefficients_mask[:, :, i] = current_correlation > \
                                                         config['badchannel_detection']['gel_bridge']['threshold']

            # Now define a new matrix with sequences of correlations above threshold
            correlation_sequence = numpy.empty((raw_data.shape[1], raw_data.shape[1]))
            for irow in range(correlation_sequence.shape[0]):
                for icol in range(correlation_sequence.shape[1]):

                    # Get the current sequence
                    current_sequence = correlation_coefficients_mask[irow, icol, :]
                    correlation_sequence[irow, icol] = numpy.sum(current_sequence) / len(current_sequence)

            # Find indexes of channels with sequences indicating gel-bridge.
            row_ind, col_ind = numpy.where(
                correlation_sequence > config['badchannel_detection']['gel_bridge']['seq_threshold'])

            # If any gel_bridged badchannel, check if they are neighbours
            if len(row_ind) > 0:

                # Load the montage to assure they are neighbours (5cm for example).
                montage = raw.get_montage()
                channels_position = montage.get_positions()
                channels_position = channels_position['ch_pos']
                channels_position = channels_position.values()
                channels_position = list(channels_position)
                channels_position = numpy.array(channels_position)

                for ichannel in range(len(row_ind)):

                    # Check the distance to assure they are neighbours (5cm for example).
                    ch_pos1 = channels_position[row_ind[ichannel], :]
                    ch_pos2 = channels_position[col_ind[ichannel], :]
                    distance = numpy.linalg.norm(ch_pos1 - ch_pos2)


                    # If the channels are close enough, they are gel-bridged
                    if distance < config['badchannel_detection']['gel_bridge']['neighbour_distance']:
                        badchannels_epoched.append(montage.ch_names[row_ind[ichannel]])
                        badchannels_epoched.append(montage.ch_names[col_ind[ichannel]])


            ###### WHOLE LENGTH

            # Reshape the data for whole length
            raw_data = numpy.transpose(raw_data,(1,2,0))
            raw_data = numpy.reshape(raw_data,(raw_data.shape[0],-1))

            # Correlation
            correlation_matrix = numpy.corrcoef(raw_data)
            correlation_matrix = numpy.triu(correlation_matrix, k=1)
            correlation_coefficients_mask = correlation_matrix > config['badchannel_detection']['gel_bridge']['threshold']

            # Find indexes of channels with sequences indicating gel-bridge.
            row_ind, col_ind = numpy.where(correlation_coefficients_mask)

            # If any gel_bridged badchannel, check if they are neighbours
            if len(row_ind) > 0:

                for ichannel in range(len(row_ind)):

                    # Check the distance to assure they are neighbours (5cm for example).
                    ch_pos1 = channels_position[row_ind[ichannel], :]
                    ch_pos2 = channels_position[col_ind[ichannel], :]
                    distance = numpy.linalg.norm(ch_pos1 - ch_pos2)

                    # If the channels are close enough, they are gel-bridged
                    if distance < config['badchannel_detection']['gel_bridge']['neighbour_distance']:
                        badchannels_whole_length.append(montage.ch_names[row_ind[ichannel]])
                        badchannels_whole_length.append(montage.ch_names[col_ind[ichannel]])


            ### PRINT
            if len(badchannels_epoched) == len(badchannels_whole_length):
                print(f"{current_sub}: {len(badchannels_epoched)} - {len(badchannels_whole_length)}")
            else:
                print(f"{current_sub}: {len(badchannels_epoched)} - {len(badchannels_whole_length)}  ********  ")



def print_channels_distance():

    databases = ['AI_Mind_database', 'LEMON_database']

    for current_database in databases:

        print(current_database)

        # Init the database
        config, files, sub, ses, task = init_database(current_database)

        # current info
        current_file = files[0]
        current_sub = sub[0]
        current_ses = ses[0]
        current_task = task[0]

        # Create the subjects following AI-Mind protocol
        bids_path = bids.create_bids_path(config, current_file, current_sub, current_ses, current_task)

        # Parameters for loading EEG  recordings
        channels_to_include = config['badchannel_detection']["channels_to_include"]
        channels_to_exclude = config['badchannel_detection']["channels_to_exclude"]

        # Load the raw EEG
        raw = aimind_mne.prepare_raw(
            config,
            bids_path,
            preload=False,
            channels_to_include=channels_to_include,
            channels_to_exclude=channels_to_exclude)

        # Load the montage to assure they are neighbours (5cm for example).
        montage = raw.get_montage()
        channels_position = montage.get_positions()
        channels_position = channels_position['ch_pos']
        channels_position = channels_position.values()
        channels_position = list(channels_position)
        channels_position = numpy.array(channels_position)

        # Channels name
        ch_names = montage.ch_names

        for ichannel in range(len(channels_position)):
            for jchannel in range(len(channels_position)):

                if ichannel == jchannel:
                    continue

                # Check the distance to assure they are neighbours (5cm for example).
                ch_pos1 = channels_position[ichannel, :]
                ch_pos2 = channels_position[jchannel, :]
                distance = numpy.linalg.norm(ch_pos1 - ch_pos2)

                if distance < 0.05:
                    print(f"{ch_names[ichannel]} - {ch_names[jchannel]}: {distance} ********")
                else:
                    print(f"{ch_names[ichannel]} - {ch_names[jchannel]}: {distance}")
