# -*- coding: utf-8 -*-
"""

Find different types of badchannels.

EEG: high_impedance_badchannels, high_variance_badchannels, power_spectrum_badchannels

MEG: flat_badchannels, high_variance_badchannels, power_spectrum_badchannel

Created on March 26 2024
@author: Federico Ramirez
"""

# Imports
import scipy
import numpy as np

import etl.io.bids as bids
import etl.tools.mnetools as aimind_mne



# Modules
def high_impedance_detection(config, bids_path):
    """

    Look for badchannels based on their impedance

    :arg
    config (dict): Configuration parameters (paths, parameters, etc)
    bids_path (dict): Path to the recording

    :returns
    List of badchannels

    """

    # Reads the channels information.
    channels = bids.read_chan (bids_path)

    # Identifies the channels with high impedance (higher than 200).
    hits = channels [ 'impedance' ] > config ['badchannel_detection'][ 'high_impedance_limit' ]
    high_impedance_badchannels = list ( channels.loc [ hits, 'name' ] )

    # Remove the excluded channels
    channels_to_exclude = config['badchannel_detection']['channels_to_exclude']
    high_impedance_badchannels = [current_channel for current_channel in high_impedance_badchannels if
                                  current_channel not in channels_to_exclude]

    return high_impedance_badchannels



def high_variance_detection(config, bids_path):
    """

    Look for badchannels based on the channel variance

    :arg
    config (dict): Configuration parameters (paths, parameters, etc)
    bids_path (dict): Path to the recording
    channels_to_include (str): Channels type to work on.

    :returns
    List of badchannels

    """

    # Read the ICA information
    sobi = bids.read_sobi(bids_path,'sobi-badchannels')

    # Parameters for loading EEG  recordings
    chan_tsv = bids.read_chan (bids_path)
    channels = list ( chan_tsv [ 'name' ] )
    freq_limits = [config['badchannel_detection']['high_variance_detection_low_freq'],
                   config['badchannel_detection']['high_variance_detection_high_freq']]
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

    # If there is EOG, remove those components
    if 'eog' in sobi.labels_.keys ():
        sobi.apply ( raw, exclude = sobi.labels_ [ 'eog' ] )

    # Select the current channels
    raw.pick(channels_to_include)

    # De-mean the channels
    raw_data = raw.get_data().copy()
    mean_per_channel = raw_data.mean(axis=1)
    raw_data_demean = raw_data - mean_per_channel[:, np.newaxis]

    # Estimate the average standard deviation of each channel
    raw_data_demean_std = raw_data_demean.std(axis=1)

    # Estimate the average standard deviation of the whole recording
    average_std = raw_data_demean_std.mean()

    # Define as badchannel any channel with a deviation >3*average_std
    hits = np.flatnonzero ( raw_data_demean_std > config['badchannel_detection']['high_variance_detection_threshold'] * average_std )
    high_variance_badchannels = [ raw.ch_names [ hit ] for hit in hits ]

    return high_variance_badchannels



def power_spectrum_detection(config,bids_path):
    """

    Look for badchannels based on anomalies in the power spectrum to detect badchannels

    :arg
    config (dict): Configuration parameters (paths, parameters, etc)
    bids_path (dict): Path to the recording
    channels_to_include (str): Channels type to work on.

    :returns
    List of badchannels

    """

    # Parameters for loading EEG  recordings
    freq_limits = [config['badchannel_detection']['pow_spectrum_detection_low_freq'],
                   config['badchannel_detection']['pow_spectrum_detection_high_freq']]
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

    # Compute the power spectrum
    win_length = 8 * raw.info['sfreq']
    freq,psd = scipy.signal.welch(
        raw.get_data(),
        raw.info['sfreq'],
        nperseg=win_length)
    psd_average = psd.mean(1)

    # Define the threshold as 3 times above the average
    threshold = config['badchannel_detection']['pow_spectrum_detection_threshold'] * psd.mean(1).mean(0)

    # Select the channels above the threshold
    hits = np.flatnonzero ( psd_average > threshold )
    power_spectrum_badchannels = [ raw.ch_names [ hit ] for hit in hits ]

    return power_spectrum_badchannels



def gel_bridge_detection(config, bids_path):
    """

    Look for badchannels based on gel bridges

    :arg
    config (dict): Configuration parameters (paths, parameters, etc)
    bids_path (dict): Path to the recording

    :returns
    List of badchannels

    """

    # Create empty list to append badchannels
    gel_bridge_badchannels = []

    # Read the ICA information
    sobi = bids.read_sobi(bids_path,'sobi-badchannels')

    # Parameters for loading EEG  recordings
    freq_limits = [config['badchannel_detection']['gel_bridge_low_freq_limit'],
                   config['badchannel_detection']['gel_bridge_high_freq_limit']]
    crop_seconds = [config['badchannel_detection']['crop_seconds']]
    channels_to_include = config['badchannel_detection']["channels_to_include"]
    channels_to_exclude = config['badchannel_detection']["channels_to_exclude"]

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
        set_annotations=False)

    # If there is EOG, remove those components
    if 'eog' in sobi.labels_.keys ():
        sobi.apply ( raw, exclude = sobi.labels_ [ 'eog' ] )

    # Get the data
    raw_data = raw.get_data()

    # Estimate the correlation matrix
    correlation_coefficients = np.corrcoef(raw_data)

    # Discard simmetrical values and diagonal = 1
    correlation_coefficients = np.triu(correlation_coefficients,k=1)

    # Find indexes of corr > threshold
    row_ind,col_ind = np.where(correlation_coefficients > config['badchannel_detection']['gel_bridge_detection_threshold'])

    # If any gel_bridged badchannel, check if they are neighbours
    if len(row_ind) > 0:

        # Load the montage to assure they are neighbours (5cm for example).
        montage = raw.get_montage()
        channels_position = montage.get_positions()
        channels_position = channels_position['ch_pos']
        channels_position = channels_position.values()
        channels_position = list(channels_position)
        channels_position = np.array(channels_position)

        for ichannel in range(len(row_ind)):

            # Check the distance to assure they are neighbours (5cm for example).
            ch_pos1 = channels_position[row_ind[ichannel],:]
            ch_pos2 = channels_position[col_ind[ichannel],:]
            distance = np.linalg.norm(ch_pos1 - ch_pos2)

            # If the channels are close enough, they are gel-bridged
            if distance < config['badchannel_detection']['neighbour_distance']:

                gel_bridge_badchannels.append(montage.ch_names[row_ind[ichannel]])
                gel_bridge_badchannels.append(montage.ch_names[col_ind[ichannel]])

    return gel_bridge_badchannels



def flat_channels_detection(config, bids_path, badchannels=None):
    """

    Look for badchannels based on low amplitudes

    :arg
    config (dict): Configuration parameters (paths, parameters, etc)
    bids_path (dict): Path to the recording

    :returns
    List of badchannels

    """

    # If there is no previous badchannels, create an empty list to append badchannels
    if badchannels is None:
        badchannels = []

    # Looks in magnetometers and gradiometers separately
    for current_sensor_type in ['mag', 'grad']:

        # Parameters to load raw MEG
        freq_limits = [config['visualization_low_freq_limit'], config['visualization_high_freq_limit']]
        crop_seconds = [config['crop_seconds']]

        # Load the raw data
        raw = aimind_mne.prepare_raw(
            config,
            bids_path,
            preload=True,
            channels_to_include=[current_sensor_type],
            channels_to_exclude=[],
            freq_limits=freq_limits,
            crop_seconds=crop_seconds,
            badchannels_to_metadata=False,
            exclude_badchannels=False,
            set_annotations=False)

        # De-mean the channels
        raw_data = raw.get_data().copy()
        mean_per_channel = raw_data.mean(axis=1)
        raw_data_demean = raw_data - mean_per_channel[:, np.newaxis]

        # Estimate the average standard deviation of each channel
        raw_data_demean_std = raw_data_demean.std(axis=1)

        # Estimate the average standard deviation of the whole recording
        average_std = raw_data_demean_std.mean()

        # Define as badchannel any channel with a deviation >3*average_std
        flat_badchannels = [raw.info['ch_names'][ichannel] for ichannel in range(len(raw.info['ch_names'])) if
                            raw_data_demean_std[ichannel] < config['badchannel_detection_flat_channel_deviation'] * average_std]

        # Append to the existing badchannels
        for new_badchannel in flat_badchannels:
            badchannels.append(new_badchannel)

    return badchannels