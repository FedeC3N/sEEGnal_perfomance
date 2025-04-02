#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Each step of the ETL called on the SRM database

Created on Mon 20/05/2024

@author: Fede
"""



# Imports
import os
import sys

import mne
import mne_bids

import aimind.meeg.io.bids as bids
import aimind.meeg.tools.mne as aimind_mne
from eeglabio.utils import export_mne_epochs


def update_bids_path(func):

    def wrapper(*args):

        # Get the original config and bids_path
        config = args[0]
        original_bids_path = args[1]

        # Update the bids_path
        bids_path = mne_bids.BIDSPath(
            subject=original_bids_path.subject,
            session=original_bids_path.session,
            task=original_bids_path.task,
            datatype=original_bids_path.datatype,
            root=original_bids_path.root,
            suffix=original_bids_path.suffix,
            extension='.vhdr')

        # Badchannel detection
        func(config,bids_path)


    return wrapper



def standardize(config,bids_path):
    """

    Convert .edf to .eeg

    """
    # Needed import
    import pybv

    # Read the original file
    raw = mne.io.read_raw_edf(bids_path.fpath)

    # Convert to .eeg
    raw_data = raw.get_data()
    sfreq = raw.info['sfreq']
    ch_names = raw.info['ch_names']
    fname_base = str(bids_path.basename)[:-4]
    folder_out = bids_path.directory
    pybv.write_brainvision(
        data=raw_data,
        sfreq=sfreq,
        ch_names=ch_names,
        fname_base=fname_base,
        folder_out=folder_out,
        overwrite=True)



@update_bids_path
def badchannel_detection(config,bids_path):
    """

    Call the ETL functions to detect badchannels in SRM database

    """

    # Add the ETL functions to the path
    sys.path.append(os.path.join('..','..','..','..','TSD','aimind.sEEGnal'))

    from aimind.etl.preprocess.badchannel_detection import estimate_badchannel_component
    import aimind.etl.tools.find_badchannels as find_badchannels

    # Estimate Independent Components
    estimate_badchannel_component(config,bids_path)

    # Create an empty list to append badchannels
    badchannels = []

    # Find channels with high variance
    high_variance_badchannels_eeg = find_badchannels.high_variance_detection(config,bids_path,'eeg')
    badchannels.extend(high_variance_badchannels_eeg)

    # Find abnormal power spectrum
    power_spectrum_badchannels_eeg = find_badchannels.power_spectrum_detection(config,bids_path,'eeg')
    badchannels.extend(power_spectrum_badchannels_eeg)

    # Find channels with gel bridge
    gel_bridge_badchannels = find_badchannels.gel_bridge_detection(config,bids_path)
    badchannels.extend(gel_bridge_badchannels)

    # Remove duplicates (if any)
    badchannels = list(set(badchannels))

    # Save the results
    bids.update_badchans(bids_path,badchannels)

    # Return the search path to the original state
    sys.path.remove(os.path.join('..','..','..','..','TSD','aimind.sEEGnal'))



@update_bids_path
def artifact_detection(config,bids_path):
    """

    Call the ETL functions to detect artifacts in SRM database

    """
    # Add the ETL functions to the path
    sys.path.append(os.path.join('..','..','..','..','TSD','aimind.sEEGnal'))

    from aimind.etl.preprocess.artifact_detection import estimate_artifact_components,muscle_detection, sensor_detection, EOG_detection

    # Estimate the SOBI components to detect artifacts
    estimate_artifact_components(config,bids_path,'sobi_artifacts')

    # Muscle and sensor artifact detection.
    # Muscle
    muscle_annotations = muscle_detection(
        config,
        bids_path,
        'eeg')

    # Save the muscle annotations for better jump detection
    _ = bids.write_annot(bids_path,muscle_annotations)

    # Sensor (jumps)
    sensor_annotations = sensor_detection(
        config,
        bids_path,
        'eeg')

    # Combine the annotations
    annotations = muscle_annotations.__add__(sensor_annotations)

    # Save the annotations in BIDS format
    _ = bids.write_annot(bids_path,annotations)

    # Estimate the SOBI components to detect artifacts
    estimate_artifact_components(config,bids_path,'sobi')

    # Muscle
    muscle_annotations = muscle_detection(
        config,
        bids_path,
        'eeg')

    # Save the muscle annotations for better jump detection
    _ = bids.write_annot(bids_path,muscle_annotations)

    # Sensor (jumps)
    sensor_annotations = sensor_detection(
        config,
        bids_path,
        'eeg')

    # EOG
    # Select the frontal channels
    frontal_channels = ['Fp1','Fpz','Fp2','AF7','AF8']
    EOG_annotations = EOG_detection(
        config,
        bids_path,
        'eeg',
        frontal_channels=frontal_channels)

    # Merge all the annotations into a MNE Annotation object
    annotations = EOG_annotations.__add__(muscle_annotations).__add__(sensor_annotations)

    # Save the annotations in BIDS format
    _ = bids.write_annot(bids_path,annotations)

    # Return the search path to the original state
    sys.path.remove(os.path.join('..','..','..','..','TSD','aimind.sEEGnal'))



@update_bids_path
def export_clean(config,bids_path):
    """

    In the public database, the clean recordings are saved in .set format (EEGLab)
    Do the same with the ETL outpu

    """

    # Load the recording
    epoch_definition = {'length':4,'overlap':0,'padding':0}
    raw = aimind_mne.prepare_raw(bids_path, preload=True,
                                 badchannels_to_metadata=True, exclude_badchannels=True,
                                 set_annotations=True, epoch=epoch_definition)

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

    # Create the output
    output_fname = os.path.abspath(os.path.join('data','SRM_database','standardized','derivatives','cleaned_epochs',
                                   'sub-' + bids_path.subject,'ses-' + bids_path.session,'eeg',bids_path.basename))[:-5]
    output_fname = output_fname + '_etl_clean.set'

    # Export
    export_mne_epochs(raw,output_fname)