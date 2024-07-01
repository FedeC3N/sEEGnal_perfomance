clear
clc
restoredefaultpath

% Paths
config.path.clean_data = '../../../../data/AI_Mind_database/cleaned';
config.path.dataset = '../../../../data/AI_Mind_database/dataset';

if ~exist(config.path.dataset), mkdir(config.path.dataset),end

% Get the total number of files
files = dir(sprintf('%s/*.set',config.path.clean_data));

% Structs for the files
dataset = struct('origin',[],'sub',[],'ses',[],'task',[],'path',[],'file',[]);
dataset_index = 1;

for ifile = 1 : numel(files)
    
    % Get the number of session for the current subject
    current_file = files(ifile);
    
    % Get the metadata
    expression = 'sub-(\w*)_ses-(\w*)_task-(\w*)_eeg_(\w*)_clean.set';
    dummy = regexp(current_file.name,expression,'tokens');

    % For the EEG experts, use 'eeg_expert'
    dataset(dataset_index).origin = dummy{1}{4};
    dataset(dataset_index).sub = dummy{1}{1};
    dataset(dataset_index).ses = dummy{1}{2};
    dataset(dataset_index).task = dummy{1}{3};
    dataset(dataset_index).file = current_file.name;
    dataset(dataset_index).path = config.path.clean_data;
    
    % Update the struct index
    dataset_index = dataset_index + 1;
    
    
end

% Save the struct
outfile = sprintf('%s/AI_Mind_dataset.mat',config.path.dataset);
save(outfile,'dataset')




