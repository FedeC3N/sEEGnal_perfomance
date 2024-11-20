# coding=utf-8
VERSION = "3.0.4"
# 3.0.4 - Minor changes in final_qa log
# 3.0.3 - Updates in gel bridge badchannel
# 3.0.2 - updates in requirements.txt
# VERSION = "3.0.1"
# Some changes have been done to aimind.meeg. Adapt ETL scripts to MEEG changes.
# VERSION = "3.0.0"
# Since badchannels and muscular artifacts tend to lead the IC estimation, there should be three ICAs:
#   - First ICA for badchannels.
#   - Second ICA for muscle and jumps.
#   - Third ICA for the other artifacts and clean signal..

# VERSION = "2.7.0" - Finish the detection artifact improvment
# VERSION = "2.6.9" - Scripts are pretty now :)
# VERSION = "2.6.8" - MEG split into mag and grad for badchannel detection - bids.py modules updated - Changes in final_qa.
# VERSION = "2.6.7" - Improve SOBI time and set time limit
# VERSION = "2.6.6" - Remove the condition of EEG has to be processed before copying scan3d
# VERSION = "2.6.5" - aimind.sdk version updated
# VERSION = "2.6.4" - final_qa inconsistency fixed
# VERSION = "2.6.3" - Path bugs fixed
# VERSION = "2.6.2" - Minor changes in badchannel detection (due to merging errors)
# VERSION = "2.6.1" - Some refactoring in bad_channels, and requirements update.
# VERSION = "2.6.0"
# Change de DW structure. Now the folder are standardized - standardized/derivatives  and curated - curated/derivatives.
#      Next step (FeatureExtraction) will read files only from curated - curated/derivatives
# Define parameters in config[''] to avoid changes in the scripts.

# VERSION = "2.5.2" Lurtis - Copied scan3d to curated
# VERSION = "2.5.1" Lurtis - Changes in config and test files
# VERSION = "2.5.0" Fede - Changes dw structure. Derivatives outputs are now stored in dw/curated.
# VERSION = "2.4.1"  Lurtis - Final quality result refactoring
# VERSION = "2.4.0"  Lurtis - Problem solved with import my_mne
# VERSION = "2.3.9"  Lurtis - Artifact detection
# VERSION = "2.3.8"  - gel_bridge_detection(bids_path, badchannels=badchannels) commented
# VERSION = "2.3.7"  Lurtis - Some refactoring in artifact detection and final quality assessment
# VERSION = "2.3.6"  Lurtis - Some refactoring in Bids.py
# VERSION = "2.3.5" Lurtis - Some refactoring in bad channel detection

# VERSION = "2.3.4"
# Lurtis - Some refactoring in standardize.py and test for avoid memory  in raweep

# VERSION = "2.3.3"
# Lurtis - Refactored metrics and metadata for adapt them to MCP

# VERSION = "2.3.2"
# Lurtis - Some changes in requirements and setup.py due issue with Numpy

# VERSION = "2.3.1"
# Fixed integration of C the module for reading CNT files.

# VERSION = "2.3.0"
# A final quality assessment has been introduced

# VERSION = "2.2.3"
# Remove fex directory. Now is a new project

# VERSION = "2.2.2"
# Minor changes of folder structures

# VERSION = "2.2.1"
# Minor changes in metrics

# VERSION = "2.2.0"
# Merge steps in ETL pipeline:

# badchannel_detection (estimate_badchannel_components + badchannel_detection)
# artifacts_detection (estimate_artifacts_components + artifact_detection)

# VERSION = "2.1.1"
# 2.1.1 dependencies due numpy conflict

# VERSION = "2.1.0"
# Include 3dscan standardization

# VERSION = "2.0.6"
# Change the connectivity folder to aimind.fex.conn

# VERSION = "2.0.5"
# Include Ricardo's connectivity

# VERSION = "2.0.4"
# Change empty lists as inputs for None

# VERSION = "2.0.3"
# MEG
# Consider as a badchannel any flat channel (or with very low values) in MEG.
# EEG
# Check impedances and mark as bad channels impedances >200
# Include gel-bridge detection for bad-channels.

# VERSION = "2. 0.2"
# Apply notch filter to EEG recordings

# VERSION = "2.0.1"
# Modifty prepare_raw to use raw.picks
# Unify prepare raw and create_dummy_RawArray

# VERSION = "2.0.0"
# MEG prcessing
    # Standarize
    # Badchannels
    # Artifacts

# VERSION = "1.1.2"
# Generate metrics for all scripts
    # standardize
    # estimate_badchannel_components
    # badchannel_detection
    # estimate_artifact_components
    # artifact_detection
# Sensor-artifact detection should use different thresholds for each channel.

# VERSION = "1.1.1"
# Merge duplicated modules inside bids.py: read_sobi and read_sobi_artifacts_derivatives.
# Change the session-code-number. Example: session code for subject 1-103 session 1 must be ses-11031Z.

# VERSION = "1.1.0"
# Merge prepare_raw and get_data.
# Change the way that prepare_raw is called in all the scripts.

# VERSION = "1.0.1"
# Merge duplicated modules inside bids.py: read_annot and read_annotations_derivatives.
# Add BIDS-derivatives output to know what file has been created/modified.

# VERSION = "1.0.0"
# First 100% functional code for EEG files.


































