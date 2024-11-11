"""

Get the metainformation of each expert and export the clean EEG
recordings to "cleaned/"

Created by Fede
12/06/2024

"""

# Imports
import os
import glob

import mne
import numpy
import scipy.io as sio

import aimind.meeg.tools.mnetools as aimind_mne
from eeglabio.utils import export_mne_epochs



# Define paths
config = {}
config['path_users'] = os.path.join('data','AI_Mind_database','EEG_experts_input')
config['path_curated'] = os.path.join('data','AI_Mind_database','curated')
config['path_cleaned'] = os.path.join('data','AI_Mind_database','cleaned')

# Some parameters
config['artifacts'] = ['eog','muscle','jump','visual']

# Read the users to process
user_txt_path = os.path.join(
    config['path_users'],
    'users.txt')
with open(user_txt_path, 'r') as f:
    all_lines = f.readlines()  # Read all lines into a list
f.close()

# Remove the return
users = [line.strip() for line in all_lines]

print()
for current_user in users:

    print(f'Currently working on user {current_user}')

    # Read the files processed by the user
    trl_path = os.path.join(
        config['path_users'],
        'trl',
        current_user,
        '*.mat')
    files = glob.glob(os.path.join(
        config['path_users'],
        'trl',
        current_user,
        '*.mat'))

    # Check if the user has any file to process
    if len(files) == 0:
        print('   No files for current user.')
        continue

    # Go through each file
    for current_file in files:

        # Load the metadata
        current_trl = sio.loadmat(current_file)

        # Create the output name to check if already exists (avoid overwrite)
        output_path = os.path.abspath(os.path.join(
            'data',
            'AI_Mind_database',
            'cleaned',
            current_trl['fileinfo'][0]['file'][0][0]))[:-4]

        output_fname = output_path + f"_{current_user}_clean.set"

        if os.path.exists(output_fname):
            print(f'   Already processed file {current_file[-36:]}')
            continue

        # Read the raw file
        header_file = glob.glob(os.path.join(
            config['path_curated'],
            'sub-' + current_trl['subject'][0],
            'ses-' + current_trl['session'][0],
            'eeg',
            '*' + current_trl['task'][0] + '*.vhdr'
        ))
        raw = mne.io.read_raw_brainvision(header_file[0],preload=True)
        raw = mne.io.RawArray(raw.get_data(),raw.info)

        # Remove some channels
        channels_to_remove = ['EMG1','EMG2','CLAV','VEOGL']
        raw.drop_channels(channels_to_remove,on_missing='ignore')

        # Set montage
        raw.set_montage('standard_1005',on_missing='ignore')

        # Set reference
        raw.set_eeg_reference()

        # Remove 50 Hz noise (and harmonics)
        raw.notch_filter(numpy.array([50,100,150,200,250]))

        # Read the badchannels
        if len(current_trl['chaninfo'][0]['bad'][0]) > 0:
            num_badchennls = len(current_trl['chaninfo'][0]['bad'][0][0])
            badchannels = [current_trl['chaninfo'][0]['bad'][0][0][i][0] for i in range(num_badchennls)]

            # Remove the badchannels
            raw.drop_channels(badchannels)

        # Read the annotations
        for current_artifact in config['artifacts']:
            try:
                start_sample = current_trl['artinfo'][0]['artifact'][0][current_artifact][0][0]['artifact'][0][0][:,0]
                end_sample = current_trl['artinfo'][0]['artifact'][0][current_artifact][0][0]['artifact'][0][0][:,1]
                duration_samples = end_sample - start_sample

                start_seconds = start_sample / raw.info['sfreq']
                duration_seconds = duration_samples / raw.info['sfreq']

                # Create the annotation class
                if 'annotations' in locals():
                    current_annotations = mne.Annotations(start_seconds,duration_seconds,current_artifact)
                    annotations = annotations.__add__(current_annotations)
                else:
                    annotations = mne.Annotations(start_seconds,duration_seconds,current_artifact)

            except:
                continue

        # Set the Annotations
        raw.set_annotations(annotations)

        # Apply clean components (brain and others)
        # Read the components and create the ICA class
        mixing = current_trl['compinfo'][0]['SOBI'][0]['EEG'][0][0][0][0]['mixing']
        unmixing = current_trl['compinfo'][0]['SOBI'][0]['EEG'][0][0][0][0]['unmixing']
        sobi = aimind_mne.build_bss(mixing,unmixing,raw.info['ch_names'])

        # Select the clean components
        components_to_include = current_trl['compinfo'][0]['SOBI'][0]['EEG'][0][0][0][0]['type']
        components_to_include = [i for i in range(len(components_to_include)) if components_to_include[i] == 0]

        # If desired components, apply them.
        if len(components_to_include) > 0:
            # Remove the eog components
            sobi.apply(raw,include=components_to_include)

        # Define epochs
        epoch = {'length':4,'overlap':0,'padding':0}

        # Create epochs
        raw = aimind_mne.get_epochs(raw,annot=raw.annotations,length=epoch['length'],overlap=epoch['overlap'],
                         padding=epoch['padding'],preload=True)

        # Export
        export_mne_epochs(raw,output_fname)
