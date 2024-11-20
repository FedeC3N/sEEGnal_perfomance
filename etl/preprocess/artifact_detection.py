# -*- coding: utf-8 -*-
"""

This module looks for different EEG/MEG artifacts.
First, it estimates the SOBI components (after excluding badchannels) and detects muscle and jumps artifacts. These
artifacts tend to affect SOBI estimation of other components.
After marking muscle and jumps artifacts, it estimates again SOBI components taking the artifacts out.
Then looks for EOG artifacts.

Created on September 27 08:43 2023
@author: Federico Ramirez

"""

# Imports
import os
import itertools
import traceback
from datetime import datetime as dt, timezone

import mne
import json
import mne_icalabel as iclabel

import aimind.etl.tools.find_artifacts as find_artifacts
import aimind.meeg.io.bids as bids
import aimind.meeg.tools.mnetools as aimind_mne



# Set the output levels
mne.utils.set_log_level(verbose='ERROR')


# Modules
def artifact_detection(config, eeg_task):
    """

    Checks if it is a EEG or MEG file and call the correspondent function for artifact detection.

    :arg
    config (dict): Configuration parameters (paths, parameters, etc)
    eeg_task (dict): Metadata of the session_id to process

    :returns
    A dict with the result of the process

    """

    if eeg_task["task_type_id"] == "eeg" or eeg_task["task_type_id"] == "meg":

        # Group Tasks by eeg_task_id EO,EC,AR, etc.
        # each eeg_task has several files (.cnt, .evt). Just use cnt.
        file_results = []
        task_data = sorted(eeg_task["eegmetadata"], key=lambda d: d['eeg_task_id'])
        for eeg_task_id, data_gen in itertools.groupby(task_data, key=lambda r: r["eeg_task_id"]):

            eeg_task_data = list(data_gen)

            # For each task, extract the current file metadata
            for file_data in eeg_task_data:

                # If it is a EEG or MEG recording (not metadata file)
                if file_data["eeg_type_id"] == "cnt" or file_data["eeg_type_id"] == "fif":

                    try:

                        # Detect artifacts for the current recording
                        bids_basename = file_artifact_detection(
                            config,
                            eeg_task['session_id'],
                            file_data )

                        # Save the results
                        file_results.append({
                            'result': 'processed',
                            'file': file_data['filename'],
                            'bids_basename': bids_basename,
                             "details": ""})

                    except Exception as e:

                        # Save the error
                        file_results.append({
                            'result': 'not_processed',
                            'file': file_data['filename'],
                            "details": f"Exception: {str(e)}, {traceback.format_exc()}"})

    else:

        raise ValueError(f"Task type {eeg_task['task_type_id']} not supported for artifact detection.")

    # Prepare the output metrics
    metrics, metadata = prepare_metrics(
        config,
        eeg_task,
        file_results)

    return {"metrics": metrics, "metadata": metadata}



def file_artifact_detection(config, session_id, file_data):
    """

    Call the artifact detection proceesses one by one for each type of recording.
    For MEG recordings, it finds artifacts in mags and grads channels separately and then merge the annotations.

    :arg
    config (dict): Configuration parameters (paths, parameters, etc)
    session_id (dict): Metadata of the session_id to process
    file_data (dict): Metadata of the file to process

    :returns
    A string with the BIDS path of the recording processed.

    """

    # First, create the BIDS path
    base_root = os.path.join(config["data_root"], config['dw_folder'], config['dw_standardized_folder'])
    bids_path = bids.build_bids (
        config,
        session_id,
        file_data,
        base_root)

    # Estimate the SOBI components to detect artifacts
    estimate_artifact_components(config, bids_path, 'sobi_artifacts')

    # Muscle and sensor artifact detection.
    # For EEG
    if bids_path.datatype == 'eeg':

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

    # For MEG
    if bids_path.datatype == 'meg':

        # Muscle
        muscle_annotations_mag = muscle_detection(
            config,
            bids_path,
            'mag')
        muscle_annotations_grad = muscle_detection(
            config,
            bids_path,
            'grad')
        muscle_annotations = muscle_annotations_mag.__add__(muscle_annotations_grad)

        # Save the muscle annotations for better jump detection
        _ = bids.write_annot(bids_path,muscle_annotations)

        # Sensor (jumps)
        sensor_annotations_mag = sensor_detection(
            config,
            bids_path,
            'mag')
        sensor_annotations_grad = sensor_detection(
            config,
            bids_path,
            'grad')
        sensor_annotations = sensor_annotations_mag.__add__(sensor_annotations_grad)

    # Combine the annotations
    annotations = muscle_annotations.__add__(sensor_annotations)

    # Save the annotations in BIDS format
    _ = bids.write_annot(bids_path, annotations)

    # Estimate the SOBI components to detect artifacts
    estimate_artifact_components(config,bids_path,'sobi')

    # I have to find again all the artifacts
    # For EEG
    if bids_path.datatype == 'eeg':

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

    # For MEG
    if bids_path.datatype == 'meg':

        # Muscle
        muscle_annotations_mag = muscle_detection(
            config,
            bids_path,
            'mag')
        muscle_annotations_grad = muscle_detection(
            config,
            bids_path,
            'grad')
        muscle_annotations = muscle_annotations_mag.__add__(muscle_annotations_grad)

        # Save the muscle annotations for better jump detection
        _ = bids.write_annot(bids_path,muscle_annotations)

        # Sensor (jumps)
        sensor_annotations_mag = sensor_detection(
            config,
            bids_path,
            'mag')
        sensor_annotations_grad = sensor_detection(
            config,
            bids_path,
            'grad')
        sensor_annotations = sensor_annotations_mag.__add__(sensor_annotations_grad)

        # EOG
        # Select the mags frontal channels
        frontal_channels = ['MEG0111','MEG0121','MEG1411','MEG1421']
        EOG_annotations_mag = EOG_detection(
            config,
            bids_path,
            'mag',
            frontal_channels=frontal_channels)
        # Select the grads frontal channels
        frontal_channels = ['MEG0112','MEG0122','MEG1412','MEG1422',
                            'MEG0113','MEG0123','MEG1413','MEG1423']
        EOG_annotations_grad = EOG_detection(
            config,
            bids_path,
            'grad',
            frontal_channels=frontal_channels)
        EOG_annotations = EOG_annotations_mag.__add__(EOG_annotations_grad)

    # Merge all the annotations into a MNE Annotation object
    annotations = EOG_annotations.__add__(muscle_annotations).__add__(sensor_annotations)

    # Save the annotations in BIDS format
    _ = bids.write_annot(bids_path,annotations)

    return bids_path.basename



def estimate_artifact_components(config,bids_path,derivatives_label):
    """

    Estimate ICA using SOBI algorithm after excluding badchannels. Then label each component.
    Save the mixing matrix, the unmixing matrix and the labels in BIDS format.

    :arg
        config (dict): Configuration parameters (paths, parameters, etc)
        bids_path (dict): Path to the recording

    """

    # Parameters for loading EEG or MEG recordings
    freq_limits = [config['component_estimation_low_freq_limit'],config['component_estimation_high_freq_limit']]
    resample_frequency = config['resampled_frequency_estimate_component']
    epoch_definition = {
        'length' : config['epoch_length'],
        'overlap' : config['epoch_overlap'],
        'padding' : config['epoch_padding']}

    # Specific EEG parameters
    if bids_path.datatype == 'eeg':
        channels_to_include = [config['eeg_artifact_detection_channels_included']]
        if config['eeg_channels_to_exclude'] == 'None':
            channels_to_exclude = None
        else:
            dummy = config['eeg_channels_to_exclude']
            channels_to_exclude = dummy.split()

    # Specific MEG parameters
    elif bids_path.datatype == 'meg':
        channels_to_include = [config['meg_artifact_detection_channels_included']]
        if config['meg_channels_to_exclude'] == 'None':
            channels_to_exclude = None
        else:
            channels_to_exclude = config['meg_channels_to_exclude']

    # Otherwise, not supported format
    else:
        raise ValueError(f"Data type {bids_path.datatype} not supported for artifact detection.")

    # Include artifacts or not (depending on the current SOBI estimation)
    if derivatives_label == 'sobi_artifacts':
        set_annotations = False
    else:
        set_annotations = True

    # Load raw EEG
    if bids_path.datatype == 'eeg':
        raw = aimind_mne.prepare_raw(
            bids_path,
            preload=True,
            channels_to_include=channels_to_include,
            channels_to_exclude=channels_to_exclude,
            resample_frequency=resample_frequency,
            freq_limits=freq_limits,
            badchannels_to_metadata=True,
            exclude_badchannels=True,
            set_annotations=set_annotations,
            epoch=epoch_definition)

    # Or load raw MEG
    if bids_path.datatype == 'meg':
        raw = aimind_mne.prepare_raw(
            bids_path,
            preload=True,
            channels_to_include=channels_to_include,
            channels_to_exclude=channels_to_exclude,
            freq_limits=freq_limits,
            badchannels_to_metadata=True,
            resample_frequency=resample_frequency,
            exclude_badchannels=True,
            set_annotations=set_annotations,
            create_dummy=True,
            epoch=epoch_definition)

    # Run SOBI
    sobi = aimind_mne.sobi(raw)

    # Label de ICs
    _ = iclabel.label_components(raw, sobi, method='iclabel')

    # Writes the SOBI data into the derivatives folder.
    _ = bids.write_sobi (bids_path, sobi, derivatives_label)



def EOG_detection(config,bids_path,channels_to_include, frontal_channels='all'):
    """

    Detect EOG artifacts

    :arg
    config (dict): Configuration parameters (paths, parameters, etc)
    bids_path (dict): Path to the recording
    channels_to_include (str): Channels type to work on.

    :returns
    MNE annotation object with detected EOG artifacts

    """

    # Find EOG aritfacts index
    EOG_index,last_sample,sfreq = find_artifacts.EOG_detection(
        config,
        bids_path,
        channels_to_include,
        frontal_channels=frontal_channels)

    # If any artifact
    if len(EOG_index) > 0:

        # Create the annotations
        EOG_annotations = find_artifacts.create_annotations(
            EOG_index,
            last_sample,
            sfreq,
            'bad_EOG')

    else:

        # Create an empty Annotation
        EOG_annotations = mne.Annotations(onset=[], duration=[], description=[])

    return EOG_annotations



def muscle_detection(config,bids_path, channels_to_include):
    """

    Detect muscle artifacts

    :arg
    config (dict): Configuration parameters (paths, parameters, etc)
    bids_path (dict): Path to the recording
    channels_to_include (str): Channels type to work on.

    :returns
    MNE annotation object with detected muscle artifacts

    """

    # Look the position of muscular artifacts
    muscle_index, last_sample, sfreq = find_artifacts.muscle_detection(
        config,
        bids_path,
        channels_to_include)

    # If any index
    if len(muscle_index) > 0:

        # Create the annotations
        muscle_index.sort()
        muscle_annotations = find_artifacts.create_annotations(
            muscle_index,
            last_sample,
            sfreq,
            'bad_muscle',
            fictional_artifact_duration=0.5)

    else:

        # If no artifacts, create empty Annotation
        muscle_annotations = mne.Annotations(onset=[], duration=[], description=[])

    return muscle_annotations



def sensor_detection(config,bids_path, channels_to_include):
    """

    Detect sensor artifacts (jumps)

    :arg
    config (dict): Configuration parameters (paths, parameters, etc)
    bids_path (dict): Path to the recording
    channels_to_include (str): Channels type to work on.

    :returns
    MNE annotation object with detected sensor artifacts

    """

    # Look the position of muscular artifacts
    sensor_index,last_sample,sfreq = find_artifacts.sensor_detection(
        config,
        bids_path,
        channels_to_include)

    # Create as Annotations
    if len(sensor_index) > 0:

        # Create the annotations
        sensor_index.sort()  # First sort the list
        sensor_annotations = find_artifacts.create_annotations(
            sensor_index,
            last_sample,
            sfreq,
            'bad_jump',
            fictional_artifact_duration=0.3)

    # If no artifacts, create empty Annotation
    else:
        sensor_annotations = mne.Annotations(onset=[], duration=[], description=[])

    return sensor_annotations



def prepare_metrics(config,eeg_task, file_results):
    """

    Prepare the metrics and results

    :arg
    eeg_task (dict): Metadata of the session_id to process
    file_results (dict): Output of the artifact detection for all recordings

    :returns
    Metrics and metadata

    """

    # Initialize the output
    general_result = "ok"
    num_files = num_processed = num_errors = 0

    # Check if each file has been processed
    for item in file_results:

        num_files += 1

        if item["result"] == "processed":
            num_processed += 1

        else:
            num_errors += 1
            general_result = "error"

    # All files proccesed: ok
    if general_result == "ok":
        message = f"Process {eeg_task['task_type_id']} artifact_detection on {eeg_task['session_id']} ended fine."

    # Else: error
    else:
        message = (f"Process {eeg_task['task_type_id']} artifact_detection on {eeg_task['session_id']} ended with errors. "
                   f"See log details.")

    # Save metrics and metadata
    metrics = {
        "result": general_result,
        "number_of_files": num_files,
        "files_processed": num_processed,
        "files_with_errors": num_errors}

    metadata = {
        "task_id": eeg_task["id"],
        "session_id": eeg_task["session_id"],
        "task_type_id": eeg_task["task_type_id"],
        "date": dt.now(timezone.utc),
        "detail": json.dumps(file_results),
        "message": message}

    return metrics, metadata
