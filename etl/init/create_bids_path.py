#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
script to create each subject

Created on Thu 17/05/2024

@author: Fede
"""


# Imports
import sys

import mne_bids



def create_bids_path(config,current_file,current_sub,current_ses,current_task):
    """

    Call the correct init function

    """

    # Build the function name
    func = f"create_{config['database']}"
    to_init = getattr(sys.modules[__name__],func)

    # Call the function
    bids_path = to_init(config,current_file,current_sub,current_ses,current_task)

    return bids_path


def create_AI_Mind_database(config,current_file,current_sub,current_ses,current_task):

    # Remove the unallowed characters
    current_sub = current_sub.replace('-','')
    current_task = current_task.replace('-','')

    # Builds the BIDS path from the metadata.
    bids_path = mne_bids.BIDSPath(
        subject=current_sub,
        session=current_ses,
        task=current_task,
        datatype='eeg',
        root=config['path']['data_root'],
        suffix='eeg',
        extension='.vhdr')

    return bids_path



