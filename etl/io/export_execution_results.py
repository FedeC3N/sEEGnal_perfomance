import os
from scipy.io import savemat,loadmat
import aimind.meeg.io.bids as bids


def export_execution_results(config,eeg_task,elapsed_time,mem_usage,process_name):

    # Get the bids_basename
    current_bids_path = bids.build_bids(config,eeg_task['session_id'],eeg_task['eegmetadata'][0],'')
    bids_basename = current_bids_path.basename[:-5]

    # Create the metadata filename
    metadata_path = os.path.join('metadata','AI_Mind_database','trl','etl')
    metadata_filename = os.path.join(metadata_path,bids_basename + '.mat')
    if not(os.path.exists(metadata_path)):
        os.mkdir(metadata_path)

    # Create a dictionary
    if os.path.exists(metadata_filename):
        result = loadmat(metadata_filename,simplify_cells=True)
    else:
        result = {'times_seconds': {},'memory_bytes':{}}

    # Save the values
    result['times_seconds'][process_name] = elapsed_time
    result['memory_bytes'][process_name] = mem_usage

    savemat(metadata_filename,result)