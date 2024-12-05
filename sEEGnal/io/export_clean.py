#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
script to obtain clean recordings from the metadata

Created on Thu 20/11/2024

@author: Fede
"""

# Imports
import os

import sEEGnal.tools.mnetools as aimind_mne
import sEEGnal.io.bids as bids
from eeglabio.utils import export_mne_epochs

def export_clean(config,bids_path):
    """

    For each recording, export it clean

    """

    if bids_path.datatype == 'eeg':

        # Create the derivatives bids_path
        out_filepath = bids.build_derivative(bids_path,'desc-sEEGnal_clean_eeg.set')

        # Check the folder
        if not os.path.exists(os.path.dirname(out_filepath)):
            os.makedirs(os.path.dirname(out_filepath))

        # Load the recording
        epoch_definition = config['artifact_detection']['epoch_definition']
        channels_to_exclude = config['component_estimation']["channels_to_exclude"]
        raw = aimind_mne.prepare_raw(
            config,
            bids_path,
            preload=True,
            badchannels_to_metadata=True,
            exclude_badchannels=True,
            channels_to_exclude=channels_to_exclude,
            set_annotations=True,
            epoch=epoch_definition)

        # Apply the clean components
        # Read the ICA information
        sobi = bids.read_sobi(bids_path, 'sobi')

        # Keep brain and other components
        components_to_include = []
        if 'brain' in sobi.labels_.keys():
            components_to_include.append(sobi.labels_['brain'])
        if 'other' in sobi.labels_.keys():
            components_to_include.append(sobi.labels_['other'])
        components_to_include = sum(components_to_include, [])

        # If desired components, apply them.
        if len(components_to_include) > 0:
            # Remove the eog components
            sobi.apply(raw, include=components_to_include)


        # Export
        export_mne_epochs(raw,str(out_filepath))