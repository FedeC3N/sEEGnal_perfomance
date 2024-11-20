# -*- coding: utf-8 -*-
"""

This module looks for badchannels in EEG/MEG recordings.
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

import etl.io.bids as bids
import etl.tools.mnetools as aimind_mne
import etl.tools.find_badchannels as find_badchannels
from etl.tools.measure_performance import measure_performance



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
    gel-bridged channels

    :arg
    config (dict): Configuration parameters (paths, parameters, etc)
    bids_path (BIDSpath): Metadata of the file to process

    :returns
    A list with the badchannels.

    """

    # Initialzies the derivatives.
    bids.init_derivatives (bids_path)

    # Estimate Independent Components
    estimate_badchannel_component(config,bids_path)

    # Create an empty list to append badchannels
    badchannels = []

    # Find channels with high impedance
    high_impedance_badchannels = find_badchannels.high_impedance_detection(config,bids_path)
    badchannels.extend(high_impedance_badchannels)

    # Find channels with high variance
    high_variance_badchannels_eeg = find_badchannels.high_variance_detection(config,bids_path)
    badchannels.extend(high_variance_badchannels_eeg)

    # Find abnormal power spectrum
    power_spectrum_badchannels_eeg = find_badchannels.power_spectrum_detection(config,bids_path)
    badchannels.extend(power_spectrum_badchannels_eeg)

    # Find channels with gel bridge
    gel_bridge_badchannels = find_badchannels.gel_bridge_detection(config,bids_path)
    badchannels.extend(gel_bridge_badchannels)

    # Remove duplicates (if any)
    badchannels = list(set(badchannels))

    # Save the results
    bids.update_badchans (bids_path, badchannels)
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
    freq_limits = [config['component_estimation']['low_freq_limit'],
                   config['component_estimation']['high_freq_limit']
                   ]
    resample_frequency = config['badchannel_detection']['resampled_frequency_estimate_component']
    epoch_definition = config['badchannel_detection']['epoch_definition']
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

    # Writes the SOBI data into the derivatives folder.
    _ = bids.write_sobi (bids_path, sobi, 'sobi-badchannels')


