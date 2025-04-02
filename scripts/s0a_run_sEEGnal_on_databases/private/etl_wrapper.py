#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Since each database is different, each step of the ETL has to be called differently

Created on Thu 17/05/2024

@author: Fede
"""



# Imports
import os
import importlib.util
import importlib.machinery



def standardize(config, bids_path):
    """

    Call the correct databases function

    """

    # Load the module
    file_to_call = f"process_{config['database']}.py"
    file_path = os.path.abspath(os.path.join('.', 'scripts', 's0a_run_sEEGnal_on_databases',
                                             'private',file_to_call))
    loader = importlib.machinery.SourceFileLoader(file_to_call,file_path)
    spec = importlib.util.spec_from_file_location(file_to_call,file_path)
    dummy_module = importlib.util.module_from_spec(spec)
    loader.exec_module(dummy_module)

    # Call the module and the specific function
    dummy_module.standardize(config,bids_path)



def badchannel_detection(config,bids_path):
    """

    Call the correct databases function

    """

    # Load the module
    file_to_call = f"process_{config['database']}.py"
    file_path = os.path.abspath(os.path.join('.', 'scripts', 's0a_run_sEEGnal_on_databases',
                                             'private', file_to_call))
    loader = importlib.machinery.SourceFileLoader(file_to_call, file_path)
    spec = importlib.util.spec_from_file_location(file_to_call, file_path)
    dummy_module = importlib.util.module_from_spec(spec)
    loader.exec_module(dummy_module)

    # Call the module and the specific function
    dummy_module.badchannel_detection(config, bids_path)



def artifact_detection(config,bids_path):
    """

    Call the correct databases function

    """

    # Load the module
    file_to_call = f"process_{config['database']}.py"
    file_path = os.path.abspath(os.path.join('.', 'scripts', 's0a_run_sEEGnal_on_databases',
                                             'private', file_to_call))
    loader = importlib.machinery.SourceFileLoader(file_to_call, file_path)
    spec = importlib.util.spec_from_file_location(file_to_call, file_path)
    dummy_module = importlib.util.module_from_spec(spec)
    loader.exec_module(dummy_module)

    # Call the module and the specific function
    dummy_module.artifact_detection(config, bids_path)



def final_qa(config,bids_path):
    """

    Call the correct databases function

    """

    # Load the module
    file_to_call = f"process_{config['database']}.py"
    file_path = os.path.abspath(os.path.join('.', 'scripts', 's0a_run_sEEGnal_on_databases',
                                             'private', file_to_call))
    loader = importlib.machinery.SourceFileLoader(file_to_call, file_path)
    spec = importlib.util.spec_from_file_location(file_to_call, file_path)
    dummy_module = importlib.util.module_from_spec(spec)
    loader.exec_module(dummy_module)

    # Call the module and the specific function
    dummy_module.final_qa(config, bids_path)


