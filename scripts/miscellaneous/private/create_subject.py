#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
script to create each subject

Created on Thu 17/05/2024

@author: Fede
"""


# Imports
import os
import sys

import mne_bids



def create_bids_path(config,current_subject_id,current_session_id,current_task):
    """

    Call the correct init function

    """

    # Build the function name
    func = f"create_{config['database']}"
    to_init = getattr(sys.modules[__name__],func)

    # Call the function
    bids_path = to_init(config,current_subject_id,current_session_id,current_task)

    return bids_path


def create_SRM_database(config,current_subject_id,current_session_id,current_task):

    # Define the metadata to create the BIDS path
    participant_id = current_subject_id
    session_id = current_session_id
    task = current_task
    datatype = 'eeg'
    base_root = os.path.join(config['data_root'],config['dw_folder'],config['dw_standardized_folder'])
    suffix = 'eeg'
    extension = '.edf'

    # Builds the BIDS path from the metadata.
    bids_path = mne_bids.BIDSPath(
        subject=participant_id,
        session=session_id,
        task=task,
        datatype=datatype,
        root=base_root,
        suffix=suffix,
        extension=extension)

    # Check if the file exists
    if not (os.path.isfile(bids_path.fpath)):
        bids_path = []

    return bids_path



def create_AI_Mind_database(config,current_subject_id,current_session_id,current_task):
    """

    Create the bids_path

    """

    # Create the correct session_id
    control_digit = "TRWAGMYFPDXBNJZSQVHLCKE"
    dummy = current_subject_id + current_session_id
    control_digit_index = int(dummy) % 23
    real_session_id = dummy + control_digit[control_digit_index]

    # Define the metadata to create the BIDS path
    participant_id = current_subject_id
    session_id = real_session_id
    task = current_task
    datatype = 'eeg'
    base_root = os.path.join(config['data_root'],config['dw_folder'],config['dw_curated_folder'])
    suffix = 'eeg'
    extension = '.vhdr'

    # Builds the BIDS path from the metadata.
    bids_path = mne_bids.BIDSPath(
        subject=participant_id,
        session=session_id,
        task=task,
        datatype=datatype,
        root=base_root,
        suffix=suffix,
        extension=extension)

    # Check if the file exists
    if not (os.path.isfile(bids_path.fpath)):
        bids_path = []

    return bids_path