clear
clc
restoredefaultpath

% Paths
config.path.clean_data = '../../../../data/SRM_database/standardized/derivatives/cleaned_epochs';
config.path.dataset = '../../../../data/SRM_database/dataset';

if ~exist(config.path.dataset), mkdir(config.path.dataset),end

% Get the total number of files
sub_folders = dir(sprintf('%s/sub*',config.path.clean_data));
sub_folders = {sub_folders.name};

% Structs for the files
dataset = struct('database',[],'sub',[],'ses',[],'task',[],'path',[],'file',[]);
dataset_index = 1;

for isub = 1 : numel(sub_folders)
    
    % Get the number of session for the current subject
    current_sub = sub_folders{isub};
    ses_folders = dir(sprintf('%s/%s/ses*',...
        config.path.clean_data,current_sub));
    ses_folders = {ses_folders.name};
    
    % Go through each ses folder to count the files
    for ises = 1 : numel (ses_folders)
        
        % Current session
        current_ses = ses_folders{ises};
        
        % Find SRM files
        dataset(dataset_index).database = 'SRM_database';
        current_file = dir(sprintf('%s/%s/%s/eeg/*desc-epochs_eeg.set',...
            config.path.clean_data,current_sub,current_ses));
        dataset(dataset_index).sub = current_sub(5:end);
        dataset(dataset_index).ses = current_ses(5:end);
        dataset(dataset_index).task = 'resteyesc';
        dataset(dataset_index).file = current_file.name;
        dataset(dataset_index).path = sprintf('%s/%s/%s/eeg',...
            config.path.clean_data,current_sub,current_ses);
        
        % Update the struct index
        dataset_index = dataset_index + 1;
        
        % Find ETL files
        dataset(dataset_index).database = 'ETL_database';
        current_file = dir(sprintf('%s/%s/%s/eeg/*etl_clean.set',...
            config.path.clean_data,current_sub,current_ses));
        dataset(dataset_index).sub = current_sub(5:end);
        dataset(dataset_index).ses = current_ses(5:end);
        dataset(dataset_index).task = 'resteyesc';
        dataset(dataset_index).file = current_file.name;
        dataset(dataset_index).path = sprintf('%s/%s/%s/eeg',...
            config.path.clean_data,current_sub,current_ses);
        
        % Update the struct index
        dataset_index = dataset_index + 1;
        
    end
    
end

% Save the struct
outfile = sprintf('%s/SRM_dataset.mat',config.path.dataset);
save(outfile,'dataset')




