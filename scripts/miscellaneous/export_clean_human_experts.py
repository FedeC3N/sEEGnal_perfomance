#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
script to run export clean on the human experts recordings

Created on Thu 20/02/2025

@author: Fede
"""

# Imports
import re
import os
import glob
import json
import numpy as np

import mne
import mne_bids
from scipy.io import loadmat
from eeglabio.utils import export_mne_epochs

import sEEGnal.tools.mnetools as mnetools

# Functions
def struct_to_dict(mat_struct):
    """Recursively converts MATLAB struct to a nested Python dictionary"""
    import numpy as np
    if isinstance(mat_struct,np.ndarray) and mat_struct.dtype.names:
        return {field:struct_to_dict(mat_struct[field]) for field in mat_struct.dtype.names}
    elif isinstance(mat_struct,np.ndarray) and mat_struct.size == 1:
        return struct_to_dict(mat_struct.item())  # Extract single element
    else:
        return mat_struct

# configs
# Read the config dictionary
with open('./config/etl_configuration.json','r') as file:
    config = json.load(file)
config['path'] = {}
config['path']['data_root'] = os.path.join('databases','AI_Mind_database')
config['path']['experts'] = os.path.join('databases','AI_Mind_database','derivatives')

# Add extra info to the dictionary
config["channels_to_include"] = [
    'Fp1', 'Fpz', 'Fp2', 'F7', 'F3', 'Fz', 'F4', 'F8', 'FC5', 'FC1', 'FC2', 'FC6',
    'M1', 'T7', 'C3', 'Cz', 'C4', 'T8', 'M2', 'CP5', 'CP1', 'CP2', 'CP6', 'P7', 'P3',
    'Pz', 'P4', 'P8', 'POz', 'O1', 'O2', 'AF7', 'AF3', 'AF4', 'AF8', 'F5', 'F1', 'F2',
    'F6', 'FC3', 'FCz', 'FC4', 'C5', 'C1', 'C2', 'C6', 'CP3', 'CP4', 'P5', 'P1', 'P2',
    'P6', 'F9', 'PO3', 'PO4', 'F10', 'FT7', 'FT8', 'TP7', 'TP8', 'PO7', 'PO8', 'FT9',
    'FT10', 'TPP9h', 'TPP10h', 'PO9', 'PO10', 'P9', 'P10', 'AFF1', 'AFz', 'AFF2', 'FFC5h',
    'FFC3h', 'FFC4h', 'FFC6h', 'FCC5h', 'FCC3h', 'FCC4h', 'FCC6h', 'CCP5h', 'CCP3h', 'CCP4h',
    'CCP6h', 'CPP5h', 'CPP3h', 'CPP4h', 'CPP6h', 'PPO1', 'PPO2', 'I1', 'Iz', 'I2', 'AFp3h', 'AFp4h',
    'AFF5h', 'AFF6h', 'FFT7h', 'FFC1h', 'FFC2h', 'FFT8h', 'FTT9h', 'FTT7h', 'FCC1h', 'FCC2h', 'FTT8h',
    'FTT10h', 'TTP7h', 'CCP1h', 'CCP2h', 'TTP8h', 'TPP7h', 'CPP1h', 'CPP2h', 'TPP8h', 'PPO9h', 'PPO5h',
    'PPO6h', 'PPO10h', 'POO9h', 'POO3h', 'POO4h', 'POO10h', 'OI1h', 'OI2h']
config["channels_to_exclude"] = ['EMG1', 'EMG2', 'CLAV', 'VEOGL']
config['artifacts_name'] = ['eog','muscle','jump','visual']

# Get the list of humnan experts
testers = glob.glob(os.path.join(config['path']['experts'],'*'))
testers = [os.path.basename(tester) for tester in testers]
testers = [tester for tester in testers if tester not in 'sEEGnal']

# Get the list of subjects
subjects = glob.glob(os.path.join(config['path']['data_root'],'sub*'))
subjects = [os.path.basename(subject) for subject in subjects]

# Get the list of files
files = []
for current_subject in subjects:
    current_file = glob.glob(os.path.join(config['path']['data_root'],current_subject,'ses-1','eeg','*eeg.eeg'))
    files.append(current_file[0])

# Extract the info
pattern = "(sub-.*)_(ses-.*)_(task-.*)_eeg.eeg"
dummies = [os.path.basename(file) for file in files]
matches = [re.findall(pattern,dummy) for dummy in dummies]
matches,files = zip(*[(match,file) for match,file in zip(matches,files) if match])
sub = [match[0][0] for match in matches]
ses = [match[0][1] for match in matches]
task = [match[0][2] for match in matches]

# For each tester
testers = ['fede']
for current_tester in testers:

    print(current_tester)


    # For each subject, do things
    for current_index in range(len(files)):

        current_sub = sub[current_index]
        current_ses = ses[current_index]
        current_task = task[current_index]

        # Create out bidspath
        outpath = os.path.join(config['path']['experts'],current_tester,'clean',sub[current_index],ses[current_index],'eeg')
        filename = f"{current_sub}_{current_ses}_{current_task}_desc-{current_tester}_clean_eeg.set"
        outfile = os.path.join(outpath,filename)

        print(f"    {filename}")

        # Check the folder
        if not os.path.exists(os.path.dirname(outfile)):
            os.makedirs(os.path.dirname(outfile))

        if os.path.exists(outfile):
            print('       Already calculated')
            continue

        # Read the metadata
        matfile = f"{current_sub}_{current_ses}_{current_task}_desc-sketch_eeg.mat"
        current_metadata = loadmat(os.path.join(outpath,matfile))

        # Remove meta fields (like '__header__', '__version__', '__globals__')
        current_metadata = {key:value for key,value in current_metadata.items() if not key.startswith("__")}

        # Process each struct in the .mat file
        current_metadata = {key:struct_to_dict(value) for key,value in current_metadata.items()}

        # Create the bidspath
        bids_path = mne_bids.BIDSPath(
            subject=current_sub[-4:],
            session=current_ses[-1:],
            task=current_task[-3:],
            datatype='eeg',
            root=config['path']['data_root'],
            suffix='eeg',
            extension='.vhdr')

        # Load the original recording
        raw = mne.io.read_raw_brainvision(bids_path,preload=True)

        # Set reference
        raw.set_eeg_reference()

        # Remove 50 Hz noise (and harmonics)
        raw.notch_filter(config['component_estimation']['notch_frequencies'])

        # Include and exclude channels explicitly
        raw.pick(config['channels_to_include'])
        raw.drop_channels(config['channels_to_exclude'],on_missing='ignore')

        # Create annotations
        annotations = []
        for current_artifact in config['artifacts_name']:

            # If no artifacts, continue
            if len(current_metadata['artinfo']['artifact'][current_artifact]['artifact']) == 0:
                continue

            onset = current_metadata['artinfo']['artifact'][current_artifact]['artifact'][:,0] / raw.info['sfreq']
            duration = (current_metadata['artinfo']['artifact'][current_artifact]['artifact'][:,1] - current_metadata['artinfo']['artifact'][current_artifact]['artifact'][:,0])/raw.info['sfreq']
            if len(annotations) > 0:
                current_annotation = mne.Annotations(onset,duration,f"bad_{current_artifact}")
                annotations = annotations.__add__(current_annotation)
            else:
                annotations = mne.Annotations(onset,duration,f"bad_{current_artifact}")


        # Add annotatons
        if len(annotations) > 0:
            raw.set_annotations(annotations)

        # Read the badchannels
        if len(current_metadata['chaninfo']['bad']) > 0:

            # Check if one badchannel or several
            if isinstance(current_metadata['chaninfo']['bad'][0],(str,)):
                badchannels = [current_metadata['chaninfo']['bad']]
            else:
                badchannels = [current_metadata['chaninfo']['bad'][0][i][0] for i in
                        range(len(current_metadata['chaninfo']['bad'][0]))]

            # Add the badchannel information
            raw.info['bads'] = badchannels

            # Exclude badchannels
            raw = raw.pick(None,exclude='bads')

        # Crop
        raw.crop(tmin=10,tmax=raw.times[-1] - 10)

        # Convert to epochs
        raw = mnetools.get_epochs(raw,annot=raw.annotations,length=4,overlap=0,padding=2,preload=True)

        # If nothing left, exclude:
        if len(raw) == 0:
            print('       No clean segments')
            continue

        # Read SOBI
        matlab_sobi = current_metadata['compinfo']['SOBI']

        # Create SOBI object
        sobi = mnetools.build_bss(matlab_sobi['mixing'],matlab_sobi['unmixing'],raw.info['ch_names'])

        # Apply SOBI
        components_to_exclude = [i for i in range(len(matlab_sobi['type'])) if matlab_sobi['type'][i] == 1]
        sobi.apply(raw,exclude=components_to_exclude)

        # Export
        export_mne_epochs(raw,outfile)

