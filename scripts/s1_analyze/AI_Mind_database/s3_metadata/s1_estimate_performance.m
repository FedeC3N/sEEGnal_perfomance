%{

Estimate the performance of sEEGnal and human experts by estimating the
time expent and memory used, and the number of rejected channels, ICs and
artifacts.


@author: Fede

%}

clear
clc
restoredefaultpath

% Add functions to the path
addpath('../shared/')

% Paths
config.path.derivatives = fullfile('..','..','..','..','databases',...
    'AI_Mind_database','derivatives');
config.path.performance = fullfile(config.path.derivatives,'sEEGnal',...
    'performance');

% Create output folder
if ~exist(config.path.performance), mkdir(config.path.performance),end

% Avoid overwrite
config.overwrite = true;

% Load the dataset
dataset_path = fullfile(config.path.derivatives,'sEEgnal','sEEGnal_dataset.mat');
load(dataset_path);

% Processes
modules = {'standardize', 'badchannel_detection', 'artifact_detection'  };
measures = {'times_seconds','memory_bytes','badchannels','artifacts','ICs_rejected'};
for ifile = 1 : numel(dataset)
    
    outfile_name = sprintf('sub-%s_ses-%s_task-%s_%s_performance.mat',...
        dataset(ifile).sub,dataset(ifile).ses,dataset(ifile).task,...
        dataset(ifile).origin);
    outfile = fullfile(config.path.performance,outfile_name);
    
    % Check if exist and overwrite
    if ~exist(outfile) || (exist(outfile) && config.overwrite)
        
        % Save the output structure
        performance = [];
        
        %%%% METADATA
        % Get the file to read
        current_file = fullfile(config.path.derivatives,'sEEGnal', 'clean',...
            ['sub-' dataset(ifile).sub],'ses-1','eeg',...
            ['*' dataset(ifile).task '*measure_performance.json']);
        current_file = dir(current_file);
        
        % Read the json file
        fid = fopen(fullfile(current_file.folder, current_file.name));
        raw = fread(fid);
        str = char(raw');
        fclose(fid);
        metadata = jsondecode(str);
        performance.times_seconds.badchannel_detection = metadata.badchannel_detection.times_seconds;
        performance.times_seconds.artifact_detection = metadata.artifact_detection.times_seconds;
        performance.memory_bytes.badchannel_detection = metadata.badchannel_detection.memory_bytes;
        performance.memory_bytes.artifact_detection = metadata.artifact_detection.memory_bytes;

        %%%% BADCHANNELS
        % Get the file to read
        current_file = fullfile(config.path.derivatives,'sEEGnal', 'clean',...
            ['sub-' dataset(ifile).sub],'ses-1','eeg',...
            ['*' dataset(ifile).task '*channels.tsv']);
        current_file = dir(current_file);
        
        % Read the file and extract the badchannels
        metadata = tdfread(fullfile(current_file.folder, current_file.name),'tab');
        badchannels = cellstr(metadata.status);
        badchannels_index = strcmp('bad',badchannels);
        badchannels = cellstr(metadata.x0xEF0xBB0xBFname);
        badchannels = badchannels(badchannels_index);
        badchannels_description = cellstr(metadata.status_description);
        badchannels_description = badchannels_description(badchannels_index);
        performance.badchannels.label = badchannels;
        performance.badchannels.description = badchannels_description;
        
        %%%% ARTIFACTS
        % Get the file to read
        current_file = fullfile(config.path.derivatives,'sEEGnal', 'clean',...
            ['sub-' dataset(ifile).sub],'ses-1','eeg',...
            ['*' dataset(ifile).task '*desc-artifacts_annotations.tsv']);
        current_file = dir(current_file);
        
        % Read the file and extract the artifacts
        try
            metadata = tdfread(fullfile(current_file.folder, current_file.name),'tab');
            performance.artifacts.label = cellstr(metadata.label);
            performance.artifacts.onset = metadata.x0xEF0xBB0xBFonset;
            performance.artifacts.duration = metadata.duration;
        catch
            performance.artifacts.label = nan;
            performance.artifacts.onset = nan;
            performance.artifacts.duration = nan;
        end
        
        %%%% ICs
        % Get the file to read
        current_file = fullfile(config.path.derivatives,'sEEGnal', 'clean',...
            ['sub-' dataset(ifile).sub],'ses-1','eeg',...
            ['*' dataset(ifile).task '*desc-sobi_artifacts_annotations.tsv']);
        current_file = dir(current_file);
        
        % Read the file and extract the ICs
        metadata = tdfread(fullfile(current_file.folder, current_file.name),'tab');
        ICs_rejected = cellstr(metadata.channel);
        ICs_rejected_labels = cellstr(metadata.label);
        performance.ICs_rejected.IC = ICs_rejected;
        performance.ICs_rejected.label = ICs_rejected_labels;   
        
        % Save the file
        save(outfile,'-struct','performance')
          
        % Add the info to the dataset
        performance = [];
        performance.path = config.path.performance;
        performance.file = outfile_name;
        dataset(ifile).performance = performance;
        
    else
        
        fprintf('   Already calculated. Do not overwrite.\n\n')
        
    end
    
end

% Save the dataset
save('-v7.3',dataset_path,'dataset')
