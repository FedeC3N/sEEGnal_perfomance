# Imports
import os
import re
import random

import aimind.meeg.io.bids as bids
import aimind.meeg.tools.mne as aimind_mne
from aimind.sdk.appconfig.client import AIMindAppConfigClient
from aimind.etl.preprocess.artifact_detection import sensor_detection


def plot_jump_artifacts(config, bids_path):

    # Detect muscle artifacts
    sensor_annotations = sensor_detection(
        config,
        bids_path,
        'eeg')

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
        resample_frequency=500,
        badchannels_to_metadata=True,
        exclude_badchannels=True,
        set_annotations=False)

    # Check if there is EOG and EKG components
    components_to_include = []
    if 'brain' in sobi.labels_:
        components_to_include.append(sobi.labels_['brain'])
    if 'other' in sobi.labels_:
        components_to_include.append(sobi.labels_['other'])

    components_to_include = sum(components_to_include,[])

    # Remove the eog components
    if len(components_to_include) > 0:
        sobi.apply(raw,include=components_to_include)

    # Add the detected annotations
    raw.set_annotations(sensor_annotations)

    # Plot
    raw.plot(block=True)



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
    session_id = '1-291-1-P'

    # Select one random recording
    subject_folder = os.path.join(folder_subjects,session_id,session_id)
    subject_files = os.listdir(subject_folder)

    # Avoid 0-AR
    not_ok = 1
    while not_ok:
        random_file_index = random.randrange(0,len(subject_files) - 1)
        current_file = subject_files[random_file_index]

        # Create the meta-information
        pattern = "[0-9]-[0-9]{3}-[0-9]-[A-Z]_(1-EO|2-EC|3-EO|4-EC)_.*"
        pattern = "[0-9]-[0-9]{3}-[0-9]-[A-Z]_([0-9]-EC)_.*"

        dummy = re.findall(pattern,current_file)
        if len(dummy) > 0:
            not_ok = 0

    eeg_task_id = dummy[0]
    file_data = {'eeg_task_id':eeg_task_id,'eeg_type_id':'cnt'}

    # Create the BIDS path
    destination_root = os.path.join(config['data_root'],config['dw_folder'],config['dw_standardized_folder'])
    bids_path = bids.build_bids(config,session_id,file_data,destination_root)

    return bids_path



# Init parameters in config
config = init_params()

# Update the config information for ETL_perfomance
config['data_root'] = '..\\..\\ETL_performance'
config['dw_folder'] = 'data'
config['datatype'] = 'eeg'

# Create the bids_path
bids_path = create_subject_bids_path(config)

# Plots
print('Working on ' + str(bids_path))

# Plot clean recordings with muscle annotations
config['channels_to_exclude'] = ['EMG1', 'EMG2', 'CLAV', 'VEOGL']
plot_jump_artifacts(config, bids_path)