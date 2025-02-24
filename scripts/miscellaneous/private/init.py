#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
script to prepare the databases to run ETL

Created on Thu 17/05/2024

@author: Fede
"""




# Imports
import os
import sys
import glob
import json

from aimind.sdk.appconfig.client import AIMindAppConfigClient


def init_config(selected_database):
    """

    Load the saved configuration file

    """

    with open('./config/etl_configuration.json','r') as file:
        config = json.load(file)
    config['path'] = {}
    config['path']['data_root'] = os.path.join('databases','AI_Mind_database')
    config['path']['experts'] = os.path.join('databases','AI_Mind_database','derivatives')
    config['database'] = selected_database
    config['tasks'] = ['1EO','2EC','3EO','4EC']
    config['artifacts_name'] = ['eog','muscle','jump','visual']

    if selected_database == 'human_experts_database':
        config['testers'] = ['fede','isa','luis','maria']
    else:
        config['testers'] = ['sEEGnal']

    return config



def init_database(config):
    """

    Call the correct databases function

    """

    # Build the function name
    func = f"init_{config['database']}"
    to_init = getattr(sys.modules[__name__],func)

    # Call the function
    subjects_id,session_id = to_init(config)

    return subjects_id,session_id



def init_SRM_database(config):

    # Folders to find the subjects
    config['data_root'] = os.getcwd()
    config['dw_folder'] = os.path.join('data',config['database'])
    folder_subjects = os.path.join(config['data_root'],config['dw_folder'],
                                   config['dw_standardized_folder'],'sub*')

    # Get all the subjects
    subjects_id = glob.glob(folder_subjects)
    subjects_id = [current[-3:] for current in subjects_id]

    # Pre-defined sessions
    sessions_id = ['t1','t2']


    return subjects_id,sessions_id



def init_AI_Mind_database(config):

    # Folders to find the subjects
    config['data_root'] = os.getcwd()
    config['dw_folder'] = os.path.join('databases',config['database'])
    folder_subjects = os.path.join(config['data_root'],config['dw_folder'],'sub*')

    # Get all the subjects
    subjects_id = glob.glob(folder_subjects)
    subjects_id = [current[-4:] for current in subjects_id]

    # Pre-defined sessions
    sessions_id = ['1','2','3','4']


    return subjects_id,sessions_id


def init_human_experts_database(config):

    # Folders to find the subjects
    config['data_root'] = os.getcwd()
    config['dw_folder'] = os.path.join('databases','AI_Mind_database')
    folder_subjects = os.path.join(config['data_root'],config['dw_folder'],'sub*')

    # Get all the subjects
    subjects_id = glob.glob(folder_subjects)
    subjects_id = [current[-4:] for current in subjects_id]

    # Pre-defined sessions
    sessions_id = ['1','2','3','4']


    return subjects_id,sessions_id