#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 20 12:40:25 2023

@author: bru
"""


import re
import mne
import numpy
import datetime

from aimind.etl.tools.old import bss as bss,my_signal as signal

'''
from aimind.etl.io.bids_old import read_badchannels_derivatives, read_annotations_derivatives
'''
import aimind.etl.io.old.bids as bids

# Function for two-pass filtering on MNE objects.
def filtfilt (mnedata, num = 1, den = 1, hilbert = False ):
    ''' Wrapper to apply two-pass filtering to MNE objects.'''
    
    
    # Checks if the data is a valid MNE object.
    if not ( isinstance ( mnedata, mne.epochs.Epochs ) or
            isinstance ( mnedata, mne.epochs.Epochs ) or
            isinstance ( mnedata, mne.io.RawArray ) ):
        
        print ( 'Unsuported data type.' )
        return None
    
    
    # Creates a copy of the input data.
    mnedata  = mnedata.copy ()
    
    # Gets the raw data matrix.
    rawdata  = mnedata.get_data ()
    
    # Filters the data.
    rawdata  = signal.filtfilt ( rawdata, num = num, den = den, hilbert = hilbert )
    
    # Updates the MNE object.
    mnedata._data = rawdata
    
    ## Creates a new MNE object with the filtered data.
    #mnedata   = mne.EpochsArray ( rawdata, data.info, events = data.events, verbose = False )
    
    ## Creates a new MNE object with the filtered data.
    #mnedata    = mne.io.RawArray ( rawdata, data.info, verbose = False )
    
    # Returns the MNE object.
    return mnedata


# Perform SOBI on MNE object
def sobi (mnedata, nlag = None, nsource = None ):
    ''' Wrapper to estimating SOBI componentes from MNE objects.'''
    
    
    # Checks if the data is a valid MNE object.
    if not ( isinstance ( mnedata, mne.epochs.Epochs ) or
            isinstance ( mnedata, mne.epochs.Epochs ) or
            isinstance ( mnedata, mne.io.RawArray ) ):
        
        print ( 'Unsuported data type.' )
        return None
    
    
    # Creates a copy of the input data.
    mnedata  = mnedata.copy ()
    
    # Gets the channel labels.
    chname   = mnedata.ch_names
    
    # Gets the raw data matrix.
    rawdata  = mnedata.get_data ()
    
    
    # Estimates the SOBI mixing matrix.
    mixing, unmixing  = bss.sobi (rawdata,nlag = nlag,nsource = nsource)
    
    
    # Builds the MNE ICA object.
    mnesobi  = build_bss ( mixing, unmixing, chname, method = 'sobi' )
    
    # Returns the MNE ICA object.
    return mnesobi


def build_raw ( info, data, montage = None ):
    
    
    # Lists the channels in the data.
    ch_label = info [ 'channels' ] [ 'label' ]
    
    
    # If no montage assumes the standard 10-05.
    if montage is None:
        montage  = mne.channels.make_standard_montage ( 'standard_1005' )
    
    
    # Identifies the EEG, EOG, ECG, and EMG channels.
    ind_eeg  = numpy.where ( numpy.in1d ( ch_label, montage.ch_names ) )
    ind_eog  = numpy.where ( [ re.search ( 'EOG', label ) != None for label in ch_label ] )
    ind_ecg  = numpy.where ( [ re.search ( 'CLAV', label ) != None for label in ch_label ] )
    ind_emg  = numpy.where ( [ re.search ( 'EMG', label ) != None for label in ch_label ] )
    
    # Marks all the channels as EEG.
    ch_types = numpy.array ( [ 'eeg' ] * len ( ch_label ) )
    
    # Sets the channel types.
    ch_types [ ind_eeg ] = 'eeg'
    ch_types [ ind_eog ] = 'eog'
    ch_types [ ind_ecg ] = 'ecg'
    ch_types [ ind_emg ] = 'emg'
    
    
    # Creates the MNE-Python information object.
    mneinfo  = mne.create_info ( 
        ch_names  = list ( info [ 'channels' ] [ 'label' ] ),
        sfreq     = info [ 'sample_rate' ],
        ch_types  = ch_types )
    
    # Adds the montage, if provided.
    if montage is not None:
        mneinfo.set_montage ( montage )
    
    
    # Creates the MNE-Python raw data object.
    mneraw   = mne.io.RawArray ( data.T, mneinfo, verbose = False )
    
    # Overwrites the default parameters.
    mneraw.set_meas_date ( info [ 'acquisition_time' ] )
    
    # Adds the calibration factor.
    mneraw._cals = numpy.ones ( len ( ch_label ) )
    
    # Marks the 'active' channels.
    mneraw._read_picks = [ numpy.arange ( len ( ch_label ) ) ]
    
    
    # Gets the information about the impedances, if any.
    if 'impedances' in info:
        
        # Takes only the first measurement.
        if len(info['impedances']) > 0:
            impmeta    = info [ 'impedances' ] [0]
            impedances = impmeta [ 'measurement' ]

            # Fills the extra information for MNE.
            for channel, value in impedances.items ():

                impedances [ channel ] = {
                    'imp':           value,
                    'imp_unit':      impmeta [ 'unit' ],
                    'imp_meas_time': datetime.datetime.fromtimestamp ( impmeta [ 'time' ] ) }


            # Adds the impedances to the MNE object.
            mneraw.impedances = impedances
    
    
    # Gets the annotations, if any.
    annotations = mne.Annotations (
        [ annot [ 'onset' ] for annot in info [ 'events' ] ],
        [ annot [ 'duration' ] for annot in info [ 'events' ] ],
        [ annot [ 'description' ] for annot in info [ 'events' ] ] )
    
    # Adds the annotations to the MNE object.
    mneraw.set_annotations ( annotations )
        
    
    # Returns the MNE Raw object.
    return mneraw


def build_bss ( mixing, unmixing, chname, icname = None, method = 'bss' ):
    
    
    # Gets the number of channels and components.
    nchannel = mixing.shape [0]
    nsource  = mixing.shape [1]
    
    # Generates the labels for the components, if no provided.
    if icname is None:
        
        # Creates the labels.
        icname   = [
            'IC%03d' % ( index + 1 )
            for index in range ( nsource ) ]
    
    
    # If the matrices are not square adds some dummy components.
    if nchannel != nsource:
        
        # Creates the dummy components.
        ndummy = nchannel - nsource
        mdummy = numpy.zeros ( [ nchannel, ndummy ] )
        
        # Creates the dummy labels.
        ldummy = [
            'DUM%03d' % ( index + 1 )
            for index in range ( ndummy ) ]
        
        # Concatenates the matrix and the labels.
        mixing   = numpy.append ( mixing, mdummy, 1 )
        unmixing = numpy.append ( unmixing, mdummy.T, 0 )
        icname   = icname + ldummy
    
    
    # Creates a dummy MNE ICA object.
    mnebss   = mne.preprocessing.ICA ()
    
    # Fills the object with the SOBI data.
    mnebss.ch_names         = chname
    mnebss._ica_names       = icname
    mnebss.mixing_matrix_   = mixing
    mnebss.unmixing_matrix_ = unmixing
    
    # Fills some dummy metadata.
    mnebss.current_fit             = 'raw'
    mnebss.method                  = method
    mnebss.n_components_           = mixing.shape[0]
    mnebss.n_iter_                 = 0
    mnebss.n_samples_              = 0
    mnebss.pca_components_         = numpy.eye ( nchannel )
    mnebss.pca_explained_variance_ = numpy.ones ( nchannel )
    mnebss.pca_mean_               = numpy.zeros ( nchannel )
    mnebss.pre_whitener_           = 0.01 * numpy.ones ( [ nchannel, 1 ] )
    
    
    # Returns the MNE ICA object.
    return mnebss


# Function to segment and MNE raw data object avoiding the artifacts.
def get_epochs(raw, annot=None, length=4, overlap=None, padding=None, preload=False):
    # Gets the annotations from the raw data, if required.
    if annot is None:
        annot = raw.annotations

    # Removes the annotations extending beyond the data.
    annot.crop(tmin=raw.times[0], tmax=raw.times[-1])

    # Sanitizes the input.
    overlap = overlap or 0
    padding = padding or 0

    # Gets the first and last valid times.
    ftime = raw.times[0] + padding
    ltime = raw.times[-1] - padding - overlap

    # Gets the separation between epochs.
    step = length - overlap

    # Gets the beginning and ending indexes for each artifact.
    artbeg = annot.onset
    artend = annot.onset + annot.duration

    # Merges overlapping artifacts.
    for index in range(len(artbeg) - 1, 0, -1):

        # Checks if the artifacts overlap.
        if artend[index - 1] > artbeg[index]:
            # Extendes the previous artifact, if required.
            artend[index - 1] = numpy.max(artend[index - 1: index + 1])

            # Removes the current artifact.
            artbeg = numpy.delete(artbeg, index)
            artend = numpy.delete(artend, index)

    # Adds virtual artifacts at the beginning and end of the data.
    artbeg = numpy.append(artbeg, ltime)
    artend = numpy.insert(artend, 0, ftime)

    # Lists the clean segments of data.
    segbegs = list(artend)
    seglens = list(artbeg - artend)

    # Initializes the list of epoch onsets.
    onsets = []

    # Iterates through each clean segment.
    for index in range(len(segbegs)):
        # Splits the segment in as many epochs as possible.
        nepoch = numpy.floor((seglens[index] - length) / step).astype(int) + 1
        sides = ((seglens[index] - length) - step * (nepoch - 1)) / 2

        # Creates the markers for the epochs.
        onset = segbegs[index] + sides + step * numpy.array(range(nepoch))

        # Adds the epochs to the list.
        onsets = onsets + list(onset)

    # Creates a series of events from the defined epochs.
    events = numpy.zeros([len(onsets), 3], 'int')
    events[:, 0] = raw.time_as_index(onsets)
    events[:, 2] = 1

    # Generates a MNE epoch structure from the data and the events.
    epochs = mne.Epochs(raw, events, tmin=-padding, tmax=length + padding, baseline=None, verbose=False,
                        preload=preload)

    # Returns the epochs object.
    return epochs


# Function to prepare MNE raw data
def prepare_raw(bids_path, preload=True, channels_to_include=None, channels_to_exclude=None,
                freq_limits=None, crop_seconds=None, badchannels_to_metadata=True, exclude_badchannels=False,
                set_annotations=True, epoch=None, create_dummy=False, resample_frequency=False):

    if channels_to_include is None:
        channels_to_include = ['all']
    if channels_to_exclude is None:
        channels_to_exclude = []

    if create_dummy:
        # Read the file
        dummy_raw = mne.io.read_raw_fif(bids_path, preload=preload)

        # Create RawArray
        dummy_raw = mne.io.RawArray(dummy_raw.get_data(), dummy_raw.info)

        # Select channels to include
        dummy_raw.pick(channels_to_include, exclude=channels_to_exclude)

        # Create the dummy data
        meg_data = dummy_raw.get_data()

        # Create the info object
        ch_names = dummy_raw.info['ch_names']
        sfreq = dummy_raw.info['sfreq']
        new_info = mne.create_info(ch_names=ch_names, sfreq=sfreq, ch_types='eeg')

        # Create the dummy RawArray
        raw = mne.io.RawArray(meg_data, new_info)

        # Create the montage for the MEG system
        ch_pos = {}
        for current_channel in dummy_raw.info['chs']:
            ch_pos[current_channel['ch_name']] = current_channel['loc'][:3]
        new_montage = mne.channels.make_dig_montage(ch_pos=ch_pos)

        # Set montage
        raw.set_montage(new_montage)

    elif bids_path.datatype == 'eeg':
        # Read the file
        raw = mne.io.read_raw_brainvision(bids_path, preload=preload)

        # Create RawArray
        raw = mne.io.RawArray(raw.get_data(), raw.info)

        # Set montage
        raw.set_montage('standard_1005', on_missing='ignore')

        # Set reference
        raw.set_eeg_reference()

        # Remove 50 Hz noise (and harmonics)
        raw.notch_filter(numpy.array([50, 100, 150, 200, 250]))

    elif bids_path.datatype == 'meg':
        # Read the file
        raw = mne.io.read_raw_fif(bids_path, preload=preload)

        # Create RawArray
        raw = mne.io.RawArray(raw.get_data(), raw.info)
    else:
        raise ValueError(f"Unsupported data type {bids_path.datatype}")

    # Include only type provided
    if not create_dummy:
        raw.pick(channels_to_include, exclude=channels_to_exclude)

    # First thing to do is epoch the data
    # if (len(epoch.keys()) == 3):
    if epoch:
        if ('length' in epoch.keys()) and ('overlap' in epoch.keys()) and ('padding' in epoch.keys()):
            raw = get_epochs(raw,annot=raw.annotations,length=epoch['length'],overlap=epoch['overlap'],
                                    padding=epoch['padding'],preload=preload)
        else:
            raise ValueError("Dictionary definition must be {'length':[],'overlap':[],'padding':[]}")

    # Downsample
    if resample_frequency:
        raw.resample(resample_frequency)

    # Filter the data
    if freq_limits:
        if len(freq_limits) != 2:
            raise ValueError("You need to define two frequency limits to filter")
            # print('You need to define two frequency limits to filter')
        else:
            raw.filter(freq_limits[0], freq_limits[1])

    # Remove the padding after filtering
    if epoch:
        raw.crop(tmin=0,tmax=raw.times[-1] - epoch['padding'])

    # Remove the beggining and the end of the recording
    if crop_seconds:
        if len(crop_seconds) == 1:
            raw.crop(tmin=crop_seconds[0], tmax=raw.times[-1] - crop_seconds[0])
        elif len(crop_seconds) == 2:
            raw.crop(tmin=crop_seconds[0], tmax=raw.times[-1] - crop_seconds[1])
        else:
            raise ValueError("The crop vector must have one or two elements")
            # print('The crop vector must have one or two elements')

    # Add the badchannel information
    if badchannels_to_metadata:

        chan = bids.read_chan (bids_path)
        badchannels = list ( chan.loc [ chan [ 'status' ] == 'bad' ] [ 'name' ] )
        '''
        badchannels = read_badchannels_derivatives(bids_path)
        '''

        # To avoid errors, the badchannels has to be among the channels included in the recording
        badchannels = list(set(badchannels) & set(raw.info['ch_names']))
        raw.info['bads'] = badchannels

        if exclude_badchannels:
            # To avoid errors if all channels are badchannels (it will be handled in final_qa).
            if not(len(badchannels) == len(raw.info['ch_names'])):
                raw = raw.pick(None, exclude='bads')

    if set_annotations:
        # Reads the annotations.
        annotations = bids.read_annot (bids_path)

        # Adds the annotations to the MNE object.
        if len ( annotations ):
            raw.set_annotations ( annotations )
        '''
        # Create the annotations from file
        onset, duration, label = read_annotations_derivatives(bids_path)
        if len(onset) > 0:
            # Add the annotations
            annotations = mne.Annotations(onset, duration, label)
            raw.set_annotations(annotations)
        '''
    return raw
