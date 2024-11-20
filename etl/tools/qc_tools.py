# -*- coding: utf-8 -*-
"""

Help functions for final check of the ETL pipeline.

Created on April 29 12:57 2024
@author: Federico Ramirez

"""


# Imports
import os
import re
import glob
import shutil
from pathlib import Path

import mne_bids

import aimind.meeg.io.bids as bids
import aimind.meeg.tools.mnetools as aimind_mne


def check_eeg_meg_derivatives_existence(config, session_id, file_data):
    """

    Check that session_id has all the expected files

    :arg
    config (dict): Configuration parameters (paths, parameters, etc)
    session_id (dict): Metadata of the session_id to process
    file_data (dict): Metadata of the file to process

    :returns
    Results of errors

    """

    # Create a list of files not found
    error_files = False
    files_not_found = []
    files_corrupted = []

    # Create the BIDS path
    base_root = os.path.join(config["data_root"], config['dw_folder'], config['dw_standardized_folder'])
    bids_path = bids.build_bids (config, session_id, file_data, base_root)

    # Standardize/
    # Find _coordsystem.json
    pattern = '*coordsystem.json'
    file_to_check = os.path.join(bids_path.directory, pattern)
    found = glob.glob(file_to_check)

    # If found, try to open the file
    if len(found) > 0:
        try:
            f = open(found[0], 'r')
            f.close()

        except Exception:
            files_corrupted.append(pattern)

    else:
        files_not_found.append(pattern)

    # Find _electrodes.tsv
    if bids_path.datatype == 'eeg':
        pattern = '*electrodes.tsv'
        file_to_check = os.path.join(bids_path.directory, pattern)
        found = glob.glob(file_to_check)

        # If found, try to open the file
        if len(found) > 0:
            try:
                f = open(found[0], 'r')
                f.close()

            except Exception:
                files_corrupted.append(pattern)

        else:
            files_not_found.append(pattern)

    # Find _channels.tsv
    pattern = '*' + bids_path.task + '_channels.tsv'
    file_to_check = os.path.join(bids_path.directory, pattern)
    found = glob.glob(file_to_check)

    # If found, try to open the file
    if len(found) > 0:
        try:
            f = open(found[0], 'r')
            f.close()

        except Exception:
            files_corrupted.append(pattern)

    else:
        files_not_found.append(pattern)

    # Decide if .eeg
    if bids_path.datatype == 'eeg':
        pattern = '*' + bids_path.task + '_eeg.eeg'

    # Decide if .fif
    elif bids_path.datatype == 'meg':
        pattern = '*' + bids_path.task + '_meg.fif'

    else:
        raise ValueError(f"Data type {bids_path.datatype} not supported for final quality assessment.")

    # Look for the file
    file_to_check = os.path.join(bids_path.directory, pattern)
    found = glob.glob(file_to_check)

    # If found, try to open the file
    if len(found) > 0:
        try:
            raw = aimind_mne.prepare_raw(bids_path, preload=False)

        except Exception:
            files_corrupted.append(found[0])

    else:
        files_not_found.append(pattern)

    # Find .json
    if bids_path.datatype == 'eeg':
        pattern = '*' + bids_path.task + '_eeg.json'

    elif bids_path.datatype == 'meg':
        pattern = '*' + bids_path.task + '_meg.json'

    else:
        raise ValueError(f"Data type {bids_path.datatype} not supported for final quality assessment.")

    # Look for the file
    file_to_check = os.path.join(bids_path.directory, pattern)
    found = glob.glob(file_to_check)

    # If found, try to open the file
    if len(found) > 0:
        try:
            f = open(found[0], 'r')
            f.close()

        except Exception:
            files_corrupted.append(pattern)

    else:
        files_not_found.append(pattern)

    # Find .vhdr
    if bids_path.datatype == 'eeg':
        pattern = '*' + bids_path.task + '_eeg.vhdr'
        file_to_check = os.path.join(bids_path.directory, pattern)
        found = glob.glob(file_to_check)

        if len(found) > 0:
            # Try to open the file
            try:
                f = open(found[0], 'r')
                f.close()
            except Exception:
                files_corrupted.append(pattern)
        else:
            files_not_found.append(pattern)

    # derivatives
    # Look for _channels.tsv
    bids_derivative_path_badchannels_tsv = bids.build_derivative(bids_path, 'channels.tsv')

    # If found, try to open the file
    if os.path.exists(bids_derivative_path_badchannels_tsv):
        try:
            _ = bids.read_chan (bids_path)

        except Exception:
            files_corrupted.append('_channels.tsv')

    else:
        files_not_found.append('_channels.tsv')

    # Look for artifacts_annotations.tsv
    bids_derivative_path_artifacts_annotations_tsv = bids.build_derivative(bids_path,
                                                                      'desc-artifacts_annotations.tsv')

    # If found, try to open the file
    if os.path.exists(bids_derivative_path_artifacts_annotations_tsv):
        try:
            _ = bids.read_annot (bids_path)

        except Exception:
            files_corrupted.append('_artifacts_annotations.tsv')

    else:
        files_not_found.append('_artifacts_annotations.tsv')

    # If any file was not found or it was corrupted, mark as error
    if len(files_corrupted) > 0:
        error_files = True
    if len(files_not_found) > 0:
        error_files = True

    return error_files,files_corrupted, files_not_found



def check_number_of_badchannels(config,session_id,file_data):
    """

    Check if the recording exceed the numbers of badchannels allowed

    :arg
    config (dict): Configuration parameters (paths, parameters, etc)
    session_id (dict): Metadata of the session_id to process
    file_data (dict): Metadata of the file to process

    :returns
    Results of errors

    """

    # Create the BIDS path
    base_root = os.path.join(config["data_root"], config['dw_folder'], config['dw_standardized_folder'])
    bids_path = bids.build_bids (config, session_id, file_data, base_root)


    # If 10% of the channels are badchannels, the recording is mark as bad
    if bids_path.datatype == 'eeg':
        threshold_bad = round(126 * config['badchannel_detection_max_number'])

    elif bids_path.datatype == 'meg':
        threshold_bad = round(306 * config['badchannel_detection_max_number'])

    else:
        raise ValueError(f"Task {bids_path.datatype} not supported for prepare_metrics badchannel estimation.")

    # Get the list of badchannels
    chan = bids.read_chan (bids_path)
    badchannels = list ( chan.loc [ chan [ 'status' ] == 'bad' ] [ 'name' ] )

    # Check if there are too many badchannels
    error_badchannels = False
    if len(badchannels) > threshold_bad:
        error_badchannels = True

    return error_badchannels,badchannels



def check_number_of_artifacts(config,session_id,file_data):
    """

    Check that there is enough clean epochs to work with.

    :arg
    config (dict): Configuration parameters (paths, parameters, etc)
    session_id (dict): Metadata of the session_id to process
    file_data (dict): Metadata of the file to process

    :returns
    Results of errors

    """

    # Create the BIDS path
    base_root = os.path.join(config["data_root"], config['dw_folder'], config['dw_standardized_folder'])
    bids_path = bids.build_bids (config, session_id, file_data, base_root)

    # Parameters to load the EEG/MEG recording
    freq_limits = [config['component_estimation_low_freq_limit'], config['component_estimation_high_freq_limit']]
    resample_frequency = config['resampled_frequency_estimate_component']
    epoch_definition = {'length': config['epoch_length'],
                        'overlap': config['epoch_overlap'],
                        'padding': config['epoch_padding']}

    # Load the raw data
    raw = aimind_mne.prepare_raw(
        bids_path,
        preload=True,
        channels_to_exclude=None,
        resample_frequency=resample_frequency,
        freq_limits=freq_limits,
        badchannels_to_metadata=False,
        exclude_badchannels=False,
        set_annotations=True,
        epoch=epoch_definition)

    # Get the number of clean epochs.
    annotations = raw.get_annotations_per_epoch()
    clean_epochs = 0
    for current_annotation in annotations:
        if not(current_annotation):
            clean_epochs += 1

    # Check if there are too many badchannels
    error_artifacts = False
    if clean_epochs < config['artifact_detection_min_epochs']:
        error_artifacts = True

    return error_artifacts, clean_epochs



def copy_standardized_files_to_curated(config,session_id,file_data):
    """

    Move the files that fulfill the quality criteria

    :arg
    config (dict): Configuration parameters (paths, parameters, etc)
    session_id (dict): Metadata of the session_id to process
    file_data (dict): Metadata of the file to process

    :returns
    Results of errors

    """

    # BIDS path to the source data (standardized)
    base_root = os.path.join(config["data_root"], config['dw_folder'], config['dw_standardized_folder'])
    bids_path_standardized = bids.build_bids (config, session_id, file_data, base_root)

    # BIDS path to the destination data (curated)
    base_root = os.path.join(config["data_root"], config['dw_folder'], config['dw_curated_folder'])
    bids_path_curated = bids.build_bids (config, session_id, file_data, base_root)

    # Create the destination folder if it does not exist
    bids_path_curated.mkdir()

    # Find the files related to the current file_data
    files_pattern = '*' + bids_path_standardized.task + '*'
    files_found = glob.glob(os.path.join(bids_path_standardized.directory, files_pattern))

    # Copy the files related to the current file_data
    for current_file in files_found:
        source = os.path.join(bids_path_standardized.directory, os.path.basename(current_file))
        destination = os.path.join(bids_path_curated.directory, os.path.basename(current_file))
        shutil.copy(source, destination)

    # Copy coordsystem.json
    files_pattern = '*coordsystem*'
    files_found = glob.glob(os.path.join(bids_path_standardized.directory, files_pattern))

    # Copy the files related to the current file_data
    for current_file in files_found:
        source = os.path.join(bids_path_standardized.directory, os.path.basename(current_file))
        destination = os.path.join(bids_path_curated.directory, os.path.basename(current_file))
        shutil.copy(source, destination)

    # EEG - Copy  electrode.tsv
    if bids_path_standardized.datatype == 'eeg':
        files_pattern = '*electrodes*'
        files_found = glob.glob(os.path.join(bids_path_standardized.directory, files_pattern))

        # Copy the files related to the current file_data
        for current_file in files_found:
            source = os.path.join(bids_path_standardized.directory, os.path.basename(current_file))
            destination = os.path.join(bids_path_curated.directory, os.path.basename(current_file))
            shutil.copy(source, destination)



def copy_derivatives_to_curated(config,session_id,file_data):
    """

    Move the derivatives of the files that fulfill the quality criteria

    :arg
    config (dict): Configuration parameters (paths, parameters, etc)
    session_id (dict): Metadata of the session_id to process
    file_data (dict): Metadata of the file to process

    :returns
    Results of errors

    """

    # BIDS path to the source data (standardized)
    base_root = os.path.join(config["data_root"], config['dw_folder'], config['dw_standardized_folder'])
    bids_path_standardized = bids.build_bids (config, session_id, file_data, base_root)

    # BIDS path to the source derivatives
    bids_path_standardized_derivatives = bids.build_derivative(bids_path_standardized,
                                                            'channels.tsv')
    bids_path_standardized_derivatives = os.path.split(bids_path_standardized_derivatives)[0]
    bids_path_standardized_derivatives = Path(bids_path_standardized_derivatives)

    # BIDS path to the destination data (curated)
    base_root = os.path.join(config["data_root"], config['dw_folder'], config['dw_curated_folder'])
    bids_path_curated = bids.build_bids (config, session_id, file_data, base_root)

    # BIDS path to the destination derivatives
    bids_path_curated_derivatives = bids.build_derivative(bids_path_curated,
                                                          'channels.tsv')
    bids_path_curated_derivatives = os.path.split(bids_path_curated_derivatives)[0]
    bids_path_curated_derivatives = Path(bids_path_curated_derivatives)

    # Create the destination folder if it does not exist
    os.makedirs(bids_path_curated_derivatives,exist_ok=True)

    # Find the files related to the current file_data
    files_pattern = '*' + bids_path_standardized.task + '*'
    files_found = glob.glob(os.path.join(bids_path_standardized_derivatives, files_pattern))

    # Copy the files related to the current file_data
    for current_file in files_found:
        source = os.path.join(bids_path_standardized_derivatives, os.path.basename(current_file))
        destination = os.path.join(bids_path_curated_derivatives, os.path.basename(current_file))
        shutil.copy(source, destination)

    # Add ETL version file
    # Path to the version file
    package_folder_path = os.path.split(os.path.split(__file__)[0])[0]

    # Read the last version line
    with open(os.path.join(package_folder_path, '_version.py'), 'r') as fd:
        version = re.search(r'^VERSION\s*=\s*[\'"]([^\'"]*)[\'"]',  # type: ignore
                            fd.read(), re.MULTILINE).group(1)
    if not version:
        raise RuntimeError('Sorry, can not find version information')

    # Create a tsv file with the version
    output_file = os.path.join(bids_path_curated_derivatives, 'ETL_version.tsv')

    # Create the file
    with open(output_file, "w") as text_file:
        text_file.write(version)



def add_extra_files(config):
    """

    Copy the files  associated to the general dataset

    :arg
    config (dict): Configuration parameters (paths, parameters, etc)

    :returns
    Results of errors

    """

    # Source and destination paths
    source_path = os.path.join(config["data_root"],config['dw_folder'],config['dw_standardized_folder'])
    destination_path = os.path.join(config["data_root"],config['dw_folder'],config['dw_curated_folder'])

    # Copy .bidsignore
    files_pattern = '.bidsignore'
    files_found = glob.glob(os.path.join(source_path,files_pattern))
    for current_file in files_found:
        source = os.path.join(source_path,os.path.basename(current_file))
        destination = os.path.join(destination_path,os.path.basename(current_file))
        shutil.copy(source,destination)

    # Copy README
    files_pattern = 'README'
    files_found = glob.glob(os.path.join(source_path,files_pattern))
    for current_file in files_found:
        source = os.path.join(source_path,os.path.basename(current_file))
        destination = os.path.join(destination_path,os.path.basename(current_file))
        shutil.copy(source,destination)



def check_scan3d_existence(config, scan3d_task):
    """

    Check that scan3d has all the expected files

    :arg
    config (dict): Configuration parameters (paths, parameters, etc)
    scan3d_task (dict): Metadata of the scan3d to process

    :returns
    Results of errors

    """

    # Create a list of files not found
    error_files = False
    files_not_found = []

    # Init IDs
    session_id = scan3d_task['session_id']
    participant_id = session_id[0:5]
    participant_id = participant_id.replace('-', '')
    session_id = session_id.replace('-', '')
    task_id = 'scan3d'
    standardized_root = os.path.join(config["data_root"], config['dw_folder'], config['dw_standardized_folder'])

    # Create the BIDS path
    bids_path = mne_bids.BIDSPath(
        subject=participant_id,
        session=session_id,
        task=task_id,
        root=standardized_root)
    standardized_scan3d_folder = os.path.join(bids_path.directory, task_id)

    # Look for the files
    scan3dmetadata = scan3d_task['scan3dmetadata']
    for current_data_type in scan3dmetadata:

        file_pattern = '*.' + current_data_type['scan3d_type_id']
        file_found = glob.glob(os.path.join(standardized_scan3d_folder, file_pattern))

        if len(file_found) == 0:
            files_not_found.append(file_pattern)

    if len(files_not_found) > 0:
        error_files = True

    return error_files, files_not_found



def copy_scan3d_to_curated(config, scan3d_task):
    """

    Move the files that fulfill the quality criteria

    :arg
    config (dict): Configuration parameters (paths, parameters, etc)
    scan3d_task (dict): Metadata of the scan3d to process

    :returns
    Results of errors

    """

    # Init IDs
    session_id = scan3d_task['session_id']
    participant_id = session_id[0:5]
    participant_id = participant_id.replace('-', '')
    session_id = session_id.replace('-', '')
    task_id = 'scan3d'
    standardized_root = os.path.join(config["data_root"], config['dw_folder'], config['dw_standardized_folder'])

    # Builds the BIDS path from the metadata.
    bids_path_standardized = mne_bids.BIDSPath(
        subject=participant_id,
        session=session_id,
        task=task_id,
        root=standardized_root)
    standardized_scan3d_folder = os.path.join(bids_path_standardized.directory, task_id)

    # Builds the BIDS path to curated.
    curated_root = os.path.join(config["data_root"], config['dw_folder'], config['dw_curated_folder'])
    bids_path_curated = mne_bids.BIDSPath(
        subject=participant_id,
        session=session_id,
        task=task_id,
        root=curated_root)

    # Create the folder if not exists
    bids_path_curated.mkdir()

    # Find the files
    scan3dmetadata = scan3d_task['scan3dmetadata']
    for current_data_type in scan3dmetadata:
        file_pattern = '*.' + current_data_type['scan3d_type_id']
        current_file = glob.glob(os.path.join(standardized_scan3d_folder, file_pattern))
        filename = os.path.basename(current_file[0])

        # Copy the file in the scan3d folder, create it if not exists
        source = os.path.join(standardized_scan3d_folder, filename)
        destination = os.path.join(bids_path_curated.directory, task_id, filename)
        os.makedirs(os.path.dirname(destination), exist_ok=True)
        shutil.copy(source, destination)

    # Mark the general process as ok
    general_result = "ok"

    return general_result