clear
clc
close all
restoredefaultpath

% Add functions to the path
addpath('../shared/')

% Paths
config.path.derivatives = fullfile('..','..','..','..','databases',...
    'LEMON_database','derivatives');

% To define later the pow matrix
complete_channel_labels = {'Fp1', 'Fp2', 'F7', 'F3', 'Fz', 'F4', 'F8', 'FC5', 'FC1', 'FC2', 'FC6', 'T7', 'C3', 'Cz', 'C4', 'T8', 'CP5', 'CP1', 'CP2', 'CP6', 'AFz', 'P7', 'P3', 'Pz', 'P4', 'P8', 'PO9', 'O1', 'Oz', 'O2', 'PO10', 'AF7', 'AF3', 'AF4', 'AF8', 'F5', 'F1', 'F2', 'F6', 'FT7', 'FC3', 'FC4', 'FT8', 'C5', 'C1', 'C2', 'C6', 'TP7', 'CP3', 'CPz', 'CP4', 'TP8', 'P5', 'P1', 'P2', 'P6', 'PO7', 'PO3', 'POz', 'PO4', 'PO8'};
config.complete_channel_labels = complete_channel_labels;

% Read the performance
performance = read_performance_dataset(config);

% Plot
plot_badchannels_in_head(config,performance);

plot_IC_distribution(performance);

plot_artifacts_distribution(performance);


%%%% FUNCTIONS
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
performance.badchannels = nan(numel(dummy.dataset),numel(config.complete_channel_labels));
performance.artifacts = [];
performance.ICs_rejected = nan(numel(dummy.dataset),7);

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
    
    % Store number of badchannels
    performance.badchannels(icurrent,:) = ismember(config.complete_channel_labels,current_performance.badchannels.label);
    
    % Store the types of artifacts
    performance.artifacts = cat(1,performance.artifacts,current_performance.artifacts.label);
    
    % Count the types of rejected ICs and store it
    for itype = 1 : numel(current_performance.ICs_rejected.label)
        
        count = numel(strsplit(current_performance.ICs_rejected.IC{itype},','));
        performance.ICs_rejected(icurrent,itype) = count;
        
    end
    
end

end


function plot_badchannels_in_head(config,performance)

% Save memory for the positions and colors
pos_elec = nan(numel(config.complete_channel_labels),2);
size_elec = 500*ones(numel(config.complete_channel_labels),1);

% Read the channels position file
lines = readlines('.private/elec1005.lay');

% Go through each line
for iline = 1 : size(lines,1)
    
    % Get the channel position in the original "complete_channel_labels"
    % for consistency
    current_line = strsplit(lines(iline),' ');
    
    if size(current_line,2) > 1
        current_channel = current_line(6);
        current_channel_index = ismember(config.complete_channel_labels,current_channel);
        
        % If present, save the position and color
        if sum(current_channel_index) == 1
            
            % Position
            pos_elec(current_channel_index,1) = current_line(2);
            pos_elec(current_channel_index,2) = current_line(3);
            
        end
        
    end
    
end

% Get the color 
color_elec = sum(performance.badchannels,1)/size(performance.badchannels,1)*100;

% Scatter the head with electrodes
fig = figure('WindowState', 'maximized');
scatter(pos_elec(:,1),pos_elec(:,2),size_elec,color_elec,'filled',...
    'MarkerEdgeColor','k')
axis equal off

% Create custom colormap (white to red)
numColors = 256; % Number of colors in the gradient
cmap = [linspace(1, 1, numColors)', linspace(1, 0, numColors)', linspace(1, 0, numColors)']; % White to red
colormap(cmap); % Set the custom colormap
colorbar

% Add channel labels
pos_labels = [pos_elec(:,1) - 0.015, pos_elec(:,2)+0.035];
text(pos_labels(:,1),pos_labels(:,2),config.complete_channel_labels)


end


function plot_IC_distribution(performance)

% IC labels
IC_labels = {'brain','muscle','eog','ecg','line_noise','ch_noise','other'};

% Count each type and estimate the percentage
count = [];
for itype = 1 : numel(IC_labels)
   
    dummy = round(sum(performance.ICs_rejected(:,itype))/100);
    dummy = itype*ones(1,dummy);
    count = cat(2,count,dummy);
    
end

% Scatter the head with electrodes
fig = figure('WindowState', 'maximized');
imagesc(count)
axis off
title('IC distribution')


end


function plot_artifacts_distribution(performance)

% IC labels
artifacst_labels = {'bad_EOG','bad_muscle','bad_jump'};

% Count each type and estimate the percentage
count = [];
for itype = 1 : numel(artifacst_labels)
   
    dummy = round(sum(ismember(performance.artifacts,artifacst_labels{itype}))/100);
    dummy = itype*ones(1,dummy);
    count = cat(2,count,dummy);
    
end

% Scatter the head with electrodes
fig = figure('WindowState', 'maximized');
imagesc(count)
axis off
title('Artifacts distribution')


end
