#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Each step of the ETL called on the AI-Mind database

Created on Mon 20/05/2024

@author: Fede
"""

# Imports
import os
import re
import sys
import time
import tracemalloc

from pprint import pprint

import aimind.meeg.io.bids as bids
import aimind.meeg.tools.mnetools as aimind_mne
from eeglabio.utils import export_mne_epochs
from scripts.s0a_run_ETL_on_databases.private.export_execution_results import export_execution_results



def create_eeg_task(func):

    def wrapper(*args):

        # Get the input
        config = args[0]
        current_session_id = args[1]

        # Create the eeg_task variable
        eeg_task = {
            'id':20921,
            'session_id':current_session_id,
            'task_type_id':'eeg',
            'eegdetails':[{
                'id':1454,
                'data_collector':'Mario Alberti',
                'head_circ_cm':51.0,
                'inion_nasion_cm':25.0,
                'cap_size':'L'}],
            'eegmetadata':[]}

        # Find all the eeg recordings
        current_path = os.path.join(config['dw_folder'],config['dw_raw_folder'],
                                    config['prospective_folder'],'eeg',current_session_id,current_session_id)
        current_subject_files = os.listdir(os.path.join(config['data_root'],current_path))

        # For each file complete the EEG variable
        for ifile in range(len(current_subject_files)):

            current_file = current_subject_files[ifile]
            type_id = current_file[-3:]

            dummy = re.findall(config['pattern'],current_file)
            if len(dummy) == 0:
                continue

            eeg_task_id = dummy[0]

            # Add the file to the metadata
            eegmetadata = {'id':97759,
                           'path':current_path,
                           'filename':current_file,
                           'eeg_task_id':eeg_task_id,
                           'eeg_type_id':type_id}

            eeg_task['eegmetadata'].append(eegmetadata)

        func(config,eeg_task)

    return wrapper



def create_curated_bids_path(func):

    def wrapper(*args):

        # Get the original config and eeg_task
        config = args[0]
        eeg_task = args[1]



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



def export_clean(config,eeg_task):
    """

    For each recording, check if it past the quality_check and then export it clean

    """

    for file_data in eeg_task['eegmetadata']:

        # Only for eeg
        if file_data["eeg_type_id"] == "cnt":

            # Curated folder
            base_root = os.path.join(config["data_root"],config['dw_folder'],config['dw_curated_folder'])

            # Create the bids_path
            current_bids_path = bids.build_bids(config,eeg_task['session_id'],file_data,base_root)

            # Check if it is in Curated folder
            if os.path.exists(current_bids_path.fpath):

                # Load the recording
                epoch_definition = {'length':4,'overlap':0,'padding':0}
                channels_to_exclude = ['CLAV','VEOGL','EMG1','EMG2']
                raw = aimind_mne.prepare_raw(current_bids_path, preload=True,
                                             badchannels_to_metadata=True, exclude_badchannels=True, channels_to_exclude=channels_to_exclude,
                                             set_annotations=True, epoch=epoch_definition)

                # Apply the clean components
                # Read the ICA information
                sobi = bids.read_sobi(current_bids_path, 'sobi')

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
                output_path = os.path.abspath(os.path.join(
                    'data',
                    config['database'],
                    'cleaned',
                    current_bids_path.basename))[:-5]

                output_fname = output_path + '_etl_clean.set'

                # Export
                export_mne_epochs(raw,output_fname)



@create_eeg_task
def standardize(config,eeg_task):
    """

    Call the ETL functions to standardize

    """

    # Add the ETL functions to the path
    sys.path.append(os.path.join('..','..','..','..','TSD','aimind.etl'))
    from aimind.etl.standardize.standardize import standardize_eeg_meg

    # Standardize
    # Trace memory and time
    start = time.time()
    tracemalloc.start()
    result = standardize_eeg_meg(config,eeg_task)
    end = time.time()
    base,peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    elapsed_time = end - start
    mem_usage = peak - base

    # Export the ETL performance results
    export_execution_results(config,eeg_task,elapsed_time,mem_usage,'standardization')

    # Some output
    print(result['metrics'])
    pprint(result['metadata']['detail'])
    print('')

    # Return the search path to the original state
    sys.path.remove(os.path.join('..','..','..','..','TSD','aimind.etl'))



@create_eeg_task
def badchannel_detection(config,eeg_task):
    """

    Call the ETL functions to standardize

    """

    # Add the ETL functions to the path
    sys.path.append(os.path.join('..','..','..','..','TSD','aimind.etl'))
    from aimind.etl.preprocess.badchannel_detection import badchannel_detection as etl_badchannel_detection

    # Badchannel detection
    # Trace memory and time
    start = time.time()
    tracemalloc.start()
    result = etl_badchannel_detection(config,eeg_task)
    end = time.time()
    base,peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    elapsed_time = end - start
    mem_usage = peak - base

    # Export the ETL performance results
    export_execution_results(config,eeg_task,elapsed_time,mem_usage,'badchannel_detection')

    # Some output
    print(result['metrics'])
    pprint(result['metadata']['detail'])
    print('')

    # Return the search path to the original state
    sys.path.remove(os.path.join('..','..','..','..','TSD','aimind.etl'))



@create_eeg_task
def artifact_detection(config,eeg_task):
    """

    Call the ETL functions to standardize

    """

    # Add the ETL functions to the path
    sys.path.append(os.path.join('..','..','..','..','TSD','aimind.etl'))
    from aimind.etl.preprocess.artifact_detection import artifact_detection as etl_artifact_detection

    # Standardize
    # Trace memory and time
    start = time.time()
    tracemalloc.start()
    result = etl_artifact_detection(config,eeg_task)
    end = time.time()
    base,peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    elapsed_time = end - start
    mem_usage = peak - base

    # Export the ETL performance results
    export_execution_results(config,eeg_task,elapsed_time,mem_usage,'artifact_detection')

    # Some output
    print(result['metrics'])
    pprint(result['metadata']['detail'])
    print('')

    # Return the search path to the original state
    sys.path.remove(os.path.join('..','..','..','..','TSD','aimind.etl'))



@create_eeg_task
def final_qa(config,eeg_task):
    """

    Call the ETL functions to standardize

    """

    # Add the ETL functions to the path
    sys.path.append(os.path.join('..','..','..','..','TSD','aimind.etl'))
    from aimind.etl.qc.final_quality_assessment import final_quality_assessment_eeg_meg

    # Standardize
    # Trace memory and time
    start = time.time()
    tracemalloc.start()
    result = final_quality_assessment_eeg_meg(config,eeg_task)
    end = time.time()
    base,peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    elapsed_time = end - start
    mem_usage = peak - base

    # Export the ETL performance results
    export_execution_results(config,eeg_task,elapsed_time,mem_usage,'final_qa')

    # Only export the valid ones
    if result['metrics']['result'] == 'ok':
        export_clean(config,eeg_task)

    # Some output
    print(result['metrics'])
    pprint(result['metadata']['detail'])
    print('')

    # Return the search path to the original state
    sys.path.remove(os.path.join('..','..','..','..','TSD','aimind.etl'))


