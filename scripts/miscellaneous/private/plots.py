#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Private function to plots

Created on Wed 22/05/2024

@author: Fede
"""



# Imports
import os

import mne
from scipy.io import loadmat

import sEEGnal.io.bids as bids
import sEEGnal.tools.mnetools as aimind_mne



# Module
def plot_raw(config,bids_path):

    # Read the raw
    raw = aimind_mne.prepare_raw(config,bids_path,
                      preload=True,
                      channels_to_include=[bids_path.datatype],
                      channels_to_exclude=['EMG1', 'EMG2', 'CLAV', 'VEOGL'],
                      freq_limits=[2, 45],
                      crop_seconds=[10],
                      badchannels_to_metadata=False,
                      exclude_badchannels=False,
                      set_annotations=False)

    raw.plot(block=True, duration=20)



def plot_clean(config,bids_path,current_tester):

    # Read the raw
    raw = mne.io.read_raw_brainvision(bids_path,preload=True)
    raw = mne.io.RawArray(raw.get_data(),raw.info)

    # Get badchannels depending on the database
    badchannels = get_badchannels(config,bids_path,current_tester)
    raw.info['bads'] = badchannels
    raw = raw.pick(None,exclude='bads')
    raw.drop_channels(['CLAV','VEOGL','EMG1','EMG2'],on_missing='ignore')

    # Add annotatons
    annotations = get_annotations(config,bids_path,raw,current_tester)
    if len(annotations) > 0:
        raw.set_annotations(annotations)

    # Get SOBI depending on the database
    sobi,components_to_include = get_sobi(config,bids_path,raw,current_tester)

    # If desired components, apply them.
    if len(components_to_include) > 0:
        # Remove the eog components
        sobi.apply(raw,include=components_to_include)

    # Filter and downsample
    raw.resample(500)
    raw.filter(2,45)
    raw.crop(10)

    # Plot
    raw.plot(block=True,duration=20)



def plot_saved_epochs(config,bids_path,current_tester):

    # Read the saved file
    outpath = os.path.join(config['path']['experts'],current_tester,'clean',f"sub-{bids_path.subject}",
                           f"ses-{bids_path.session}",'eeg')
    saved_file = f"sub-{bids_path.subject}_ses-{bids_path.session}_task-{bids_path.task}_desc-{current_tester}_clean_eeg.set"
    raw = mne.io.read_epochs_eeglab(os.path.join(outpath,saved_file))

    # Filter and resample
    raw.resample(500)
    raw.filter(2,45)

    # Plot
    raw.plot(block=True,n_epochs=2)


def get_badchannels(config,bids_path,current_tester):

    badchannels = []

    if config['database'] == 'human_experts_database':

        # Read the metadata
        outpath = os.path.join(config['path']['experts'],current_tester,'clean',f"sub-{bids_path.subject}",
                               f"ses-{bids_path.session}",'eeg')
        matfile = f"sub-{bids_path.subject}_ses-{bids_path.session}_task-{bids_path.task}_desc-sketch_eeg.mat"
        current_metadata = loadmat(os.path.join(outpath,matfile))

        # Remove meta fields (like '__header__', '__version__', '__globals__')
        current_metadata = {key:value for key,value in current_metadata.items() if not key.startswith("__")}

        # Process each struct in the .mat file
        current_metadata = {key:struct_to_dict(value) for key,value in current_metadata.items()}

        # Read the badchannels
        if len(current_metadata['chaninfo']['bad']) > 0:

            # Check if one badchannel or several
            if isinstance(current_metadata['chaninfo']['bad'][0],(str,)):
                badchannels = [current_metadata['chaninfo']['bad']]
            else:
                badchannels = [current_metadata['chaninfo']['bad'][0][i][0] for i in
                               range(len(current_metadata['chaninfo']['bad'][0]))]

    else:

        # Read the txt file saved
        chan = bids.read_chan(bids_path)
        badchannels = list(chan.loc[chan['status'] == 'bad']['name'])


    return badchannels



def get_annotations(config,bids_path,raw,current_tester):

    annotations = []

    if config['database'] == 'human_experts_database':

        # Read the metadata
        outpath = os.path.join(config['path']['experts'],current_tester,'clean',f"sub-{bids_path.subject}",
                               f"ses-{bids_path.session}",'eeg')
        matfile = f"sub-{bids_path.subject}_ses-{bids_path.session}_task-{bids_path.task}_desc-sketch_eeg.mat"
        current_metadata = loadmat(os.path.join(outpath,matfile))

        # Remove meta fields (like '__header__', '__version__', '__globals__')
        current_metadata = {key:value for key,value in current_metadata.items() if not key.startswith("__")}

        # Process each struct in the .mat file
        current_metadata = {key:struct_to_dict(value) for key,value in current_metadata.items()}

        for current_artifact in config['artifacts_name']:

            # If no artifacts, continue
            if len(current_metadata['artinfo']['artifact'][current_artifact]['artifact']) == 0:
                continue

            onset = current_metadata['artinfo']['artifact'][current_artifact]['artifact'][:,0] / raw.info['sfreq']
            duration = (current_metadata['artinfo']['artifact'][current_artifact]['artifact'][:,1] -
                        current_metadata['artinfo']['artifact'][current_artifact]['artifact'][:,0]) / raw.info['sfreq']
            if len(annotations) > 0:
                current_annotation = mne.Annotations(onset,duration,f"bad_{current_artifact}")
                annotations = annotations.__add__(current_annotation)
            else:
                annotations = mne.Annotations(onset,duration,f"bad_{current_artifact}")

    else:

        # Reads the annotations.
        annotations = bids.read_annot(bids_path)

    return annotations



def get_sobi(config,bids_path,raw,current_tester):

    sobi = []
    components_to_include = []

    if config['database'] == 'human_experts_database':

        # Read the metadata
        outpath = os.path.join(config['path']['experts'],current_tester,'clean',f"sub-{bids_path.subject}",f"ses-{bids_path.session}",'eeg')
        matfile = f"sub-{bids_path.subject}_ses-{bids_path.session}_task-{bids_path.task}_desc-sketch_eeg.mat"
        current_metadata = loadmat(os.path.join(outpath,matfile))

        # Remove meta fields (like '__header__', '__version__', '__globals__')
        current_metadata = {key:value for key,value in current_metadata.items() if not key.startswith("__")}

        # Process each struct in the .mat file
        current_metadata = {key:struct_to_dict(value) for key,value in current_metadata.items()}

        # Read SOBI
        matlab_sobi = current_metadata['compinfo']['SOBI']

        # Create SOBI object
        sobi = aimind_mne.build_bss(matlab_sobi['mixing'],matlab_sobi['unmixing'],raw.info['ch_names'])

        # Apply SOBI
        components_to_include = [i for i in range(len(matlab_sobi['type'])) if matlab_sobi['type'][i] == 0]

    else:

        # Read the ICA information
        sobi = bids.read_sobi(bids_path,'sobi')

        # Keep brain and other components
        components_to_include = []
        if 'brain' in sobi.labels_.keys():
            components_to_include.append(sobi.labels_['brain'])
        if 'other' in sobi.labels_.keys():
            components_to_include.append(sobi.labels_['other'])
        components_to_include = sum(components_to_include,[])

    return sobi, components_to_include



def struct_to_dict(mat_struct):
    """Recursively converts MATLAB struct to a nested Python dictionary"""
    import numpy as np
    if isinstance(mat_struct,np.ndarray) and mat_struct.dtype.names:
        return {field:struct_to_dict(mat_struct[field]) for field in mat_struct.dtype.names}
    elif isinstance(mat_struct,np.ndarray) and mat_struct.size == 1:
        return struct_to_dict(mat_struct.item())  # Extract single element
    else:
        return mat_struct