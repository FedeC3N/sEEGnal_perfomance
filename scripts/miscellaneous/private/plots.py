#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Private function to plots

Created on Wed 22/05/2024

@author: Fede
"""



# Imports
import aimind.meeg.io.bids as bids
import aimind.meeg.tools.mne as aimind_mne



# Module
def plot_raw(bids_path):

    # Read the raw
    raw = aimind_mne.prepare_raw(bids_path,
                      preload=True,
                      channels_to_include=[bids_path.datatype],
                      channels_to_exclude=['EMG1', 'EMG2', 'CLAV', 'VEOGL'],
                      freq_limits=[2, 45],
                      crop_seconds=[10],
                      badchannels_to_metadata=False,
                      exclude_badchannels=False,
                      set_annotations=False)

    raw.plot(block=True, duration=20)



def plot_clean(bids_path):

    # Read the raw
    raw = aimind_mne.prepare_raw(
        bids_path,
        preload=True,
        channels_to_include=[bids_path.datatype],
        channels_to_exclude=['CLAV','VEOGL','EMG1','EMG2'],
        freq_limits=[2,45],
        crop_seconds=[10],
        badchannels_to_metadata=True,
        resample_frequency=500,
        exclude_badchannels=False,
        set_annotations=True)


    # Read the ICA information
    sobi = bids.read_sobi(bids_path,'sobi')

    # Keep brain and other components
    components_to_include = []
    if 'brain' in sobi.labels_.keys():
        components_to_include.append(sobi.labels_['brain'])
    if 'other' in sobi.labels_.keys():
        components_to_include.append(sobi.labels_['other'])
    components_to_include = sum(components_to_include,[])

    # If desired components, apply them.
    if len(components_to_include) > 0:
        # Remove the eog components
        sobi.apply(raw,include=components_to_include)

    # Plot
    raw.plot(block=True,duration=20)
