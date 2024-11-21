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

from aimind.sdk.appconfig.client import AIMindAppConfigClient


def init_lurtis():
    """

    Load Lurtis configuration from pypi server

    """

    # Lurtis configuration
    os.environ['AI_MIND_CONNECTION_ENDPOINT'] = 'https://pms.lurtis.com'
    os.environ['AI_MIND_API_CREDENTIALS'] = '{"username": "etluser","password": "80Mind80@"}'
    SUBSYSTEM = "ETL"
    APP_NAME = "etl-meeg-test"
    app_config = AIMindAppConfigClient(APP_NAME,SUBSYSTEM,logging=False)
    config = app_config.get_current_configuration()
    logger = app_config.get_logger()
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
    config['dw_folder'] = os.path.join('data','SRM_database')
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
    config['dw_folder'] = os.path.join('data','AI_Mind_database')
    folder_subjects = os.path.join(config['data_root'],config['dw_folder'],
                                   config['dw_raw_folder'],'prospective','eeg')


    # Filter the subjects of interest
    subjects_id = os.listdir(folder_subjects)

    # Pre-defined sessions
    sessions_id = ['']


    return subjects_id,sessions_id
