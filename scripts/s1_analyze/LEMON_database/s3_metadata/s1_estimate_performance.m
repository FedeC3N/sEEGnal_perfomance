clear
clc
restoredefaultpath

% Add functions to the path
addpath('../shared/')

% Paths
config.path.derivatives = fullfile('..','..','..','..','databases',...
    'LEMON_database','derivatives');

% Avoid overwrite
config.overwrite = false;

% Load the dataset
dataset_path = fullfile(config.path.derivatives,'sEEgnal','sEEGnal_dataset.mat');
load(dataset_path);

% Processes
processes = {'standardize', 'badchannel_detection', 'artifact_detection'  };
measures = {'times_seconds','memory_bytes'};
all_data = nan(numel(dataset),numel(measures),numel(processes));
for ifile = 1 : numel(dataset)
    
    % Get the file to read
    current_file = fullfile(config.path.derivatives,'sEEGnal', 'clean',...
        ['sub-' dataset(ifile).sub],'ses-1','eeg',...
        ['*' dataset(ifile).task '*measure_performance.json']);
    current_file = dir(current_file);
    
    % Add the info to the dataset
    performance = [];
    performance.path = fullfile('databases','LEMON_database',...
        'derivatives''sEEgnal','clean',['sub-' dataset(ifile).sub],...
        'ses-1','eeg');
    performance.file = current_file.name;
    dataset(ifile).performance = performance;
    
    % Read the json file
    fid = fopen(fullfile(current_file.folder, current_file.name));
    raw = fread(fid);
    str = char(raw');
    fclose(fid);
    data = jsondecode(str);
    
    % Save the data
    for iprocess = 1 : numel(processes)
        for imeasure = 1 : numel(measures)
            
            all_data(ifile,imeasure,iprocess) = data.(processes{iprocess}).(measures{imeasure});
            
        end
    end
    
end

% Save the dataset
save('-v7.3',dataset_path,'dataset')

% Display the results
data_average = squeeze(mean(all_data,1));
data_average(2,:) = data_average(2,:)*0.000001; % Convert to MB
data_std = squeeze(std(all_data,1));
data_std(2,:) = data_std(2,:)*0.000001; % Convert to MB
for imeasure = 1 : numel(measures)
    
    fprintf('%s\n',measures{imeasure})
    
    for iprocess = 1 : numel(processes)
        
        fprintf('  %s: %.2f +- %.2f\n', ...
            processes{iprocess},...
            data_average(imeasure,iprocess),data_std(imeasure,iprocess))
        
    end
    fprintf('\n')
end
