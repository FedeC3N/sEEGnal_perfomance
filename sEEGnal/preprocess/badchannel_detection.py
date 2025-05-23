# -*- coding: utf-8 -*-
"""

This module looks for badchannels in EEG recordings.
First, it estimates the SOBI components and labels them.
Then looks for badchannels based on: amplitude of the signal, power spectrum, impedances.

Created on April 27 15:32 2023
@author: Federico Ramirez

"""

# Imports
import traceback
from datetime import datetime as dt, timezone

import mne
import mne_icalabel as iclabel

import sEEGnal.io.bids as bids
import sEEGnal.tools.mnetools as aimind_mne
import sEEGnal.tools.find_badchannels as find_badchannels
from sEEGnal.tools.find_badchannels import impossible_amplitude_detection
from sEEGnal.tools.measure_performance import measure_performance



# Set the output levels
mne.utils.set_log_level(verbose='ERROR')



# Modules
@measure_performance
def badchannel_detection(config, bids_path):
    """

    Checks if it is a EEG file and call the correspondent function for badchannel detection.

    :arg
    config (dict): Configuration parameters (paths, parameters, etc)
    bids_path (BIDSpath): Metadata to process

    :returns
    A dict with the result of the process

    """

    # For EEG
    if bids_path.datatype == 'eeg':

        try:

            # Detect badchannels in the current recording
            badchannels = eeg_badchannel_detection(config, bids_path)

            # Save the results
            now = dt.now(timezone.utc)
            formatted_now = now.strftime("%d-%m-%Y %H:%M:%S")
            results = {'result':'ok',
                       'bids_basename':bids_path.basename,
                       "date":formatted_now,
                       'badchannels': badchannels
                       }

        except Exception as e:

            # Save the error
            now = dt.now(timezone.utc)
            formatted_now = now.strftime("%d-%m-%Y %H:%M:%S")
            results = {'result':'error',
                   'bids_basename':bids_path.basename,
                   "date":formatted_now,
                   "details": f"Exception: {str(e)}, {traceback.format_exc()}"
                   }

    else:

        # Not accepted type to process
        now = dt.now(timezone.utc)
        formatted_now = now.strftime("%d-%m-%Y %H:%M:%S")
        results = {'result':'error',
                   'file':bids_path.basename,
                   'details':'Not accepted type of file to process',
                   'date':formatted_now
                   }


    return results



def eeg_badchannel_detection(config,bids_path):
    """

    Call the badchannel detection processes one by one for each type of recording.
    It will define badchannels based on the impedance, high amplitude, flat channels, high-energy pow spectrum, and
    gel-bridged channels.

    The order of detection is important:
        1. High impedance
        2. Flat channels
        3. Impossible high amplitude
        4. High power spectrum power
        5. Gel-bridged channels
        6. High variance channels

    :arg
    config (dict): Configuration parameters (paths, parameters, etc)
    bids_path (BIDSpath): Metadata of the file to process

    :returns
    A list with the badchannels.

    """

    # Initialzies the derivatives.
    bids.init_derivatives (bids_path)

    # Estimate Independent Components
    #estimate_badchannel_component(config,bids_path)

    # Create an empty list to append badchannels
    badchannels = []
    badchannels_description = []

    # Find channels with high impedance
    """high_impedance_badchannels = find_badchannels.high_impedance_detection(config, bids_path)
    badchannels.extend(high_impedance_badchannels)
    current_badchannel_description = ['bad_high_impedance' for i in range(len(high_impedance_badchannels))]
    badchannels_description.extend(current_badchannel_description)

    # Find channels with biologically impossible amplitude
    impossible_amplitude_badchannels = find_badchannels.impossible_amplitude_detection(config, bids_path,badchannels)
    badchannels.extend(impossible_amplitude_badchannels)
    current_badchannel_description = ['bad_impossible_amplitude_badchannels' for i in range(len(impossible_amplitude_badchannels))]
    badchannels_description.extend(current_badchannel_description)

    # Find abnormal power spectrum
    power_spectrum_badchannels = find_badchannels.power_spectrum_detection(config, bids_path,badchannels)
    badchannels.extend(power_spectrum_badchannels)
    current_badchannel_description = ['bad_power_spectrum' for i in range(len(power_spectrum_badchannels))]
    badchannels_description.extend(current_badchannel_description)

    # Find channels with gel bridge
    gel_bridge_badchannels = find_badchannels.gel_bridge_detection(config, bids_path,badchannels)
    badchannels.extend(gel_bridge_badchannels)
    current_badchannel_description = ['bad_gel_bridge' for i in range(len(gel_bridge_badchannels))]
    badchannels_description.extend(current_badchannel_description)"""

    # Find channels with high variance
    high_deviation_badchannels = find_badchannels.high_deviation_detection(config, bids_path, badchannels)
    badchannels.extend(high_deviation_badchannels)
    current_badchannel_description = ['bad_high_deviation' for i in range(len(high_deviation_badchannels))]
    badchannels_description.extend(current_badchannel_description)
    print(badchannels)

    # Save the results
    bids.update_badchans (bids_path, badchannels, badchannels_description)
    return badchannels



def estimate_badchannel_component(config, bids_path):
    """

    Estimate ICA using SOBI algorithm. Then label each component.
    Save the mixing matrix, the unmixing matrix and the labels in BIDS format.

    :arg
        config (dict): Configuration parameters (paths, parameters, etc)
        bids_path (dict): Path to the recording

    """

    # Parameters for loading EEG recordings
    freq_limits = [config['component_estimation']['low_freq'],
                   config['component_estimation']['high_freq']
                   ]
    resample_frequency = config['component_estimation']['resampled_frequency']
    epoch_definition = config['component_estimation']['epoch_definition']
    channels_to_include = config['component_estimation']["channels_to_include"]
    channels_to_exclude = config['component_estimation']["channels_to_exclude"]

    # Load raw EEG
    raw = aimind_mne.prepare_raw(
        config,
        bids_path,
        preload=True,
        channels_to_include=channels_to_include,
        channels_to_exclude=channels_to_exclude,
        resample_frequency=resample_frequency,
        freq_limits=freq_limits,
        badchannels_to_metadata=False,
        exclude_badchannels=False,
        set_annotations=False,
        epoch=epoch_definition)

    # Run SOBI
    sobi = aimind_mne.sobi(raw)

    # Label de ICs
    _ = iclabel.label_components(raw, sobi, 'iclabel')

    # Check the probabilities max probabilities of each category and find those under 0.7
    max_probability = sobi.labels_scores_.max(axis=1)
    index_unclear = [i for i in range(len(max_probability)) if max_probability[i] < config['component_estimation']['unclear_threshold']]

    # Re-arrange the label matrix and assign "other" to the unclears
    labels = ['brain', 'muscle', 'eog', 'ecg', 'line_noise', 'ch_noise']
    for current_label in labels:
        sobi.labels_[current_label] = [i for i in sobi.labels_[current_label] if i not in index_unclear]

    dummy = sobi.labels_['other'].copy()
    dummy.extend(index_unclear)
    dummy = list(set(dummy))
    sobi.labels_['other'] = dummy.copy()

    # Writes the SOBI data into the derivatives folder.
    _ = bids.write_sobi (bids_path, sobi, 'sobi-badchannels')


