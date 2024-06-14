# Imports
import os
import re
import random

import aimind.meeg.io.bids as bids
import aimind.meeg.tools.mne as aimind_mne
from aimind.sdk.appconfig.client import AIMindAppConfigClient
from aimind.etl.preprocess.artifact_detection import EOG_detection


def plot_VEOGL_and_clean_frontal_channels(config,bids_path):

    # Detect EOG artifacts
    frontal_channels = config['frontal_channels']
    EOG_annotations = EOG_detection(
        config,
        bids_path,
        config['channels_to_include'],
        frontal_channels=frontal_channels)

    # Read the ICA information
    sobi = bids.read_sobi(bids_path,'sobi')

    # EEG params
    freq_limits = [config['visualization_low_freq_limit'],config['visualization_high_freq_limit']]
    crop_seconds = [config['crop_seconds']]
    channels_to_exclude = config['channels_to_exclude']

    # Load EEG
    raw = aimind_mne.prepare_raw(
        bids_path,
        preload=True,
        freq_limits=freq_limits,
        crop_seconds=crop_seconds,
        channels_to_exclude=channels_to_exclude,
        badchannels_to_metadata=True,
        exclude_badchannels=True,
        set_annotations=False)

    # Check if there is EOG and EKG components
    components_to_exclude = []
    if 'eog' in sobi.labels_:
        components_to_exclude.append(sobi.labels_['eog'])
    if 'ecg' in sobi.labels_:
        components_to_exclude.append(sobi.labels_['ecg'])

    components_to_exclude = sum(components_to_exclude,[])

    # Remove the eog components
    if len(components_to_exclude) > 0:
        sobi.apply(raw,exclude=components_to_exclude)

    # Add the detected annotations
    raw.set_annotations(EOG_annotations)

    # Select frontal channels and VEOGL
    channels = sobi.ch_names
    frontal_channels = [current_channel for current_channel in frontal_channels if current_channel in channels]
    if config['datatype'] == 'eeg':
        channels_to_include = frontal_channels + ['VEOGL']
    else:
        channels_to_include = frontal_channels
    raw.pick(channels_to_include)

    raw.plot(block=True,proj=False)



def init_params():
    # Params
    os.environ['AI_MIND_CONNECTION_ENDPOINT'] = 'https://pms.lurtis.com'
    os.environ['AI_MIND_API_CREDENTIALS'] = '{"username": "etluser","password": "80Mind80@"}'
    SUBSYSTEM = "ETL"
    APP_NAME = "etl-meeg-test"
    app_config = AIMindAppConfigClient(APP_NAME,SUBSYSTEM,logging=False)
    config = app_config.get_current_configuration()

    return config



def create_subject_bids_path(config):

    # Select a session_id
    folder_subjects = os.path.join(config['data_root'],config['dw_folder'],config['dw_raw_folder'],
                                   config['prospective_folder'],config['datatype'])
    sessions_id = os.listdir(folder_subjects)
    # session_id = sessions_id[0]
    randon_subject_index = random.randrange(0,len(sessions_id) - 1)
    session_id = sessions_id[randon_subject_index]
    # session_id = '2-111-1-C'

    # Select one random recording
    subject_folder = os.path.join(folder_subjects,session_id,session_id)
    subject_files = os.listdir(subject_folder)

    # Avoid 0-AR
    not_ok = 1
    while not_ok:
        random_file_index = random.randrange(0,len(subject_files) - 1)
        current_file = subject_files[random_file_index]

        # Create the meta-information
        pattern = "[0-9]-[0-9]{3}-[0-9]-[A-Z]_(1-EO|3-EO)_.*"

        dummy = re.findall(pattern,current_file)
        if len(dummy) > 0:
            not_ok = 0

    eeg_task_id = dummy[0]
    if config['datatype'] == 'eeg':
        eeg_type_id = 'cnt'
    if config['datatype'] == 'meg':
        eeg_type_id = 'fif'
    file_data = {'eeg_task_id':eeg_task_id,'eeg_type_id':eeg_type_id}

    # Create the BIDS path
    destination_root = os.path.join(config['data_root'],config['dw_folder'],config['dw_standardized_folder'])
    bids_path = bids.build_bids(config,session_id,file_data,destination_root)

    return bids_path



# Init parameters in config
config = init_params()

# Update the config information for ETL_perfomance
#config['data_root'] = '..\\..\\ETL_performance'
#config['dw_folder'] = 'data'
config['datatype'] = 'meg'
config['channels_to_include'] = 'mag'

# Create the bids_path
bids_path = create_subject_bids_path(config)

# Plots
print('Working on ' + str(bids_path))

# VEOG vs frontal channels
# Select EEG frontal channels
if config['datatype'] == 'eeg':
    config['frontal_channels'] = ['Fp1', 'Fpz', 'Fp2', 'AF7', 'AF8']
    config['channels_to_exclude'] = ['EMG1', 'EMG2', 'CLAV']

# Select MEG frontal channels
if config['datatype'] == 'meg':
    if config['channels_to_include'] == 'mag':
        config['frontal_channels'] = ['MEG0111','MEG0121','MEG1411','MEG1421']
        config['channels_to_exclude'] = []

    if config['channels_to_include'] == 'grad':
        config['frontal_channels'] = ['MEG0112','MEG0122','MEG1412','MEG1422',
                                      'MEG0113','MEG0123','MEG1413','MEG1423']
        config['channels_to_exclude'] = []

plot_VEOGL_and_clean_frontal_channels(config,bids_path)