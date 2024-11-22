# -*- coding: utf-8 -*-
"""
Decide how to read a original EEG file based on the file type

Created on Thu 21/11/2024

@author: Fede
"""

# Imports
import mne
import etl.io.eep as eep

def read_source_files(config,source_filepath):


    # Detect the EEG file type based on the extension
    extension = source_filepath.split('.')[-1]

    # Brain Vision
    if extension == 'vhdr':

        # Read the file
        mnedata = mne.io.read_raw(source_filepath,preload=True)

        # Create RawArray
        mnedata = mne.io.RawArray(mnedata.get_data(),mnedata.info)

    elif extension == 'cnt':

        # Read the data using Ricardo's function
        mnedata = eep.read_mne(source_filepath)


    return mnedata