#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
script to run each step of the ETL proccess across all subjects (simulating Airflow call)

Created on Thu 19/09/2023

@author: Fede
"""

import os
import re
from aimind.sdk.appconfig.client import AIMindAppConfigClient
from aimind.etl.io.old.bids import create_bids_path
from aimind.etl.tools.old.my_mne_tools import prepare_raw
import random
from aimind.etl.preprocess.badchannel_detection import power_spectrum_detection

# Para importar desde el terminal
# py.exe 'C:\Users\feder\AppData\Local\Programs\Python\Python39\Scripts\import_config.py' https://pms.lurtis.com 'C:\Users\feder\OneDrive - Universidad Complutense de Madrid (UCM)\Mi Unidad\Proyectos - Programación\AI-Mind\TSD\aimind.etl\config\etl_configuration.yaml' ETL etluser '80Mind80@'


os.environ['AI_MIND_CONNECTION_ENDPOINT'] = 'https://pms.lurtis.com'
os.environ['AI_MIND_API_CREDENTIALS'] = '{"username": "etluser","password": "80Mind80@"}'

SUBSYSTEM = "ETL"
APP_NAME = "etl-meeg-test"

app_config = AIMindAppConfigClient(APP_NAME, SUBSYSTEM, logging=False)
config = app_config.get_current_configuration()
print(config)
logger = app_config.get_logger()


# Folders to find the subjects
task_type_id = 'eeg'
config['pow_spectrum_detection_threshold'] = 3
config['data_root'] = '..\\..\\ETL_performance'
config['dw_folder'] = 'data'
folder_subjects = os.path.join(config['data_root'], config['dw_folder'], config['dw_raw_folder'], config['prospective_folder'], task_type_id)

sessions_id = os.listdir(folder_subjects)
randon_subject_index = random.randrange(0,len(sessions_id)-1)
session_id = sessions_id[randon_subject_index]

# Select one random recording
subject_folder = os.path.join(folder_subjects,session_id,session_id)
subject_files = os.listdir(subject_folder)

# Avoid 0-AR
not_ok = 1
while not_ok:
    random_file_index = random.randrange(0,len(subject_files)-1)
    current_file = subject_files[random_file_index]

    # Create the meta-information
    pattern = "[0-9]-[0-9]{3}-[0-9]-[A-Z]_(1-EO|2-EC|3-EO|4-EC)_.*"

    dummy = re.findall(pattern,current_file)
    if len(dummy) > 0:
        not_ok = 0

eeg_task_id = dummy[0]
file_data = {'eeg_task_id':eeg_task_id,'eeg_type_id':'cnt'}

# Create the BIDS path
destination_root = os.path.join(config['data_root'], config['dw_folder'], config['dw_standardized_folder'])
bids_path = create_bids_path(config, session_id, file_data,destination_root)

power_spectrum_badchannels = power_spectrum_detection(config,bids_path)
print(power_spectrum_badchannels)


visualization_raw = prepare_raw(bids_path,preload=True, channels_to_exclude=['EMG1','EMG2','CLAV','VEOGL'],
                  freq_limits=[2,45],crop_seconds=[10],badchannels_to_metadata=False,
                  exclude_badchannels=False,
                  set_annotations=False,epoch={})
visualization_raw.plot()




# Plot the power spectrum
"""fig,(ax1,ax2) = plt.subplots(2,1)
ax1.plot(range(psd.shape[0]),psd_average,'o')
ax1.set_xticks(range(psd.shape[0]),channels,rotation='vertical')
ax1.set_xticklabels(channels)
channels_to_plot = channels
index_to_plot = [channels.index(i) for i in channels_to_plot]
for ichannel in index_to_plot:
    ax2.plot(freq,psd[ichannel,:],linewidth=1)
channels_to_plot = power_spectrum_badchannels
index_to_plot = [channels.index(i) for i in channels_to_plot]
for ichannel in index_to_plot:
    ax2.plot(freq,psd[ichannel,:],linewidth=2,color='red')
ax2.set_xlim([40,60])
plt.show()"""
