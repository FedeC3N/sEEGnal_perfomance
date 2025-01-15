clear
clc
restoredefaultpath

% Add functions to the path
addpath('../shared/')

% Paths
config.path.derivatives = fullfile('..','..','..','..','databases',...
    'LEMON_database','derivatives');

% Read the performance
performance = read_performance_dataset(config);

% Analyze
modules = {'standardize', 'badchannel_detection', 'artifact_detection'  };
measures = {'times_seconds','memory_bytes','badchannels','artifacts','ICs_rejected'};
for imeasure = 1 : numel(measures)
    
    % Current measure
    current_measure = measures{imeasure};
    fprintf('%s\n',current_measure);
    
    % Print
    if imeasure < 3
        
        for imodule = 1 : numel(modules)
            current_module = modules{imodule};
            data_average = mean(performance.(current_measure).(current_module));
            data_std = std(performance.(current_measure).(current_module));
            fprintf('  %s: %.2f +- %.2f\n', ...
                modules{imodule}, data_average,data_std)
        end
    else
        data_average = mean(performance.(current_measure));
        data_std = std(performance.(current_measure));
        fprintf('  Number: %.2f +- %.2f\n', ...
            data_average,data_std)
    end

    
end

% Read all the performance files
function performance = read_performance_dataset(config)

% Load the dataset
dataset_path = fullfile(config.path.derivatives,'sEEgnal','sEEGnal_dataset.mat');
dummy = load(dataset_path);

% Output structure
modules = {'standardize', 'badchannel_detection', 'artifact_detection'  };
performance = [];
for imodule = 1 : numel(modules)
    performance.times_seconds .(modules{imodule})= nan(numel(dummy.dataset),1);
    performance.memory_bytes.(modules{imodule}) = nan(numel(dummy.dataset),1);
end
performance.badchannels = nan(numel(dummy.dataset),1);
performance.artifacts = nan(numel(dummy.dataset),1);
performance.ICs_rejected = nan(numel(dummy.dataset),1);

for icurrent = 1 : numel(dummy.dataset)
    
    % Load performance
    current_dataset = dummy.dataset(icurrent);
    performance_file = fullfile(current_dataset.performance.path,...
        current_dataset.performance.file);
    current_performance = load(performance_file);
    
    % Store times and bytes
    for imodule = 1 : numel(modules)
        performance.times_seconds.(modules{imodule})(icurrent) = current_performance.times_seconds.(modules{imodule});
        performance.memory_bytes.(modules{imodule})(icurrent) = current_performance.memory_bytes.(modules{imodule}) * 0.000001;
    end
    
    % Store number of badchannels and artifacts
    performance.badchannels(icurrent) = numel(current_performance.badchannels.label);
    performance.artifacts(icurrent) = numel(current_performance.artifacts.label);
    
    % Count the number of rejected ICs and store it
    dummy_count = strjoin(current_performance.ICs_rejected.IC,',');
    dummy_count = strsplit(dummy_count,',');
    performance.ICs_rejected(icurrent) = numel(dummy_count);
    
end

end


