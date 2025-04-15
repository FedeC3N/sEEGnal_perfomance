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

% Avoid overwrite
config.overwrite = true;

% Artifacts
artifact_type = {'eog','muscle','jump','visual'};

% Get the different testers (except sEEGnal)
testers = dir(sprintf('%s/*',config.path.derivatives));
testers = testers(3:end-1);
testers = {testers.name};

for itester = 1 : numel(testers)
    
    current_tester = testers{itester};
    config.path.performance = fullfile(config.path.derivatives,current_tester,...
        'performance');
    
    % Create output folder
    if ~exist(config.path.performance), mkdir(config.path.performance),end
    
    % Load the dataset
    dataset_name = sprintf('%s_dataset.mat',current_tester);
    dataset_path = fullfile(config.path.derivatives,current_tester,dataset_name);
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
            current_file = fullfile(config.path.derivatives,current_tester, 'clean',...
                ['sub-' dataset(ifile).sub],'ses-1','eeg',...
                ['*' dataset(ifile).task '*trl*']);
            current_file = dir(current_file);
            metadata_trl = load(sprintf('%s/%s',current_file.folder,current_file.name));
            
            % Times and memory
            performance.times_seconds.badchannel_detection = metadata_trl.times_seconds.select_badchannels;
            performance.times_seconds.artifact_detection = nansum([metadata_trl.times_seconds.select_artifacts,...
                metadata_trl.times_seconds.component_extraction,metadata_trl.times_seconds.component_revision,...
                metadata_trl.times_seconds.artifact_revision]);
            performance.memory_bytes.badchannel_detection = metadata_trl.memory_bytes.select_badchannels;
            performance.memory_bytes.artifact_detection = max([metadata_trl.times_seconds.select_artifacts,...
                metadata_trl.times_seconds.component_extraction, metadata_trl.times_seconds.component_revision,...
                metadata_trl.times_seconds.artifact_revision]);
            
            %%%% BADCHANNELS
            % Read the file and extract the badchannels
            performance.badchannels.label = metadata_trl.chaninfo.bad;
            
            %%%% ARTIFACTS
            performance.artifacts.label = {};
            performance.artifacts.onset = [];
            performance.artifacts.duration = [];
            for itype = 1 : numel(artifact_type)
                current_artifact = artifact_type{itype};
                n_artifacts = size(metadata_trl.artinfo.artifact.(current_artifact).artifact,1);
                if n_artifacts > 0
                    [dummy{1:n_artifacts}] = deal(current_artifact);
                    performance.artifacts.label = cat(1,performance.artifacts.label, dummy');
                    clear dummy
                    dummy = metadata_trl.artinfo.artifact.(current_artifact).artifact(:,1)/2000;
                    performance.artifacts.onset = cat(1,performance.artifacts.onset, dummy);
                    clear dummy
                    dummy = (metadata_trl.artinfo.artifact.(current_artifact).artifact(:,2) - metadata_trl.artinfo.artifact.(current_artifact).artifact(:,1))/2000;
                    performance.artifacts.duration = cat(1,performance.artifacts.duration, dummy);
                    clear dummy
                end
            end
            
            %%%% ICs
            IC_rejected_index = find(metadata_trl.compinfo.SOBI.EEG.type > 0);
            performance.ICs_rejected.IC = IC_rejected_index;
            
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
    
    
end