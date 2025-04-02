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



def create_bids_path(config,current_subject_id,current_session_id):
    """

    Call the correct databases function

    """

    # Build the function name
    func = f"create_{config['database']}"
    to_init = getattr(sys.modules[__name__],func)

    # Call the function
    bids_path = to_init(config,current_subject_id,current_session_id)

    return bids_path


def create_SRM_database(config,current_subject_id,current_session_id):

    # Define the metadata to create the BIDS path
    participant_id = current_subject_id
    session_id = current_session_id
    task = 'resteyesc'
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



def create_AI_Mind_database(config,current_subject_id,current_session_id):
    """

    In this case we are creating a fake bids_path and then create the eegtask needed with a wrapper

    """

    # Dummy bids_path
    bids_path = current_subject_id

    return bids_path