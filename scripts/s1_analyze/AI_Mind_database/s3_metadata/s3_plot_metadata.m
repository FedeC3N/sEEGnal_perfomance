%{

Different plots to investigate the results.


@author: Fede

%}

clear
clc
close all
restoredefaultpath

% Add functions to the path
addpath('../shared/')

% Paths
config.path.derivatives = fullfile('..','..','..','..','databases',...
    'AI_Mind_database','derivatives');
config.path.figures = '../../../../docs/manuscript/figures/AI_Mind_database/metadata_results';

if ~exist(config.path.figures), mkdir(config.path.figures),end

% To define later the pow matrix
complete_channel_labels = {'Fp1', 'Fpz', 'Fp2', 'F7', 'F3', 'Fz', 'F4', 'F8', 'FC5', 'FC1', 'FC2', 'FC6',...
    'T7', 'C3', 'Cz', 'C4', 'T8', 'CP5', 'CP1', 'CP2', 'CP6', 'P7', 'P3',...
    'Pz', 'P4', 'P8', 'POz', 'O1', 'O2', 'AF7', 'AF3', 'AF4', 'AF8', 'F5', 'F1', 'F2',...
    'F6', 'FC3', 'FCz', 'FC4', 'C5', 'C1', 'C2', 'C6', 'CP3', 'CP4', 'P5', 'P1', 'P2',...
    'P6', 'F9', 'PO3', 'PO4', 'F10', 'FT7', 'FT8', 'TP7', 'TP8', 'PO7', 'PO8', 'FT9',...
    'FT10', 'TPP9h', 'TPP10h', 'PO9', 'PO10', 'P9', 'P10', 'AFF1', 'AFz', 'AFF2', 'FFC5h',...
    'FFC3h', 'FFC4h', 'FFC6h', 'FCC5h', 'FCC3h', 'FCC4h', 'FCC6h', 'CCP5h', 'CCP3h', 'CCP4h',...
    'CCP6h', 'CPP5h', 'CPP3h', 'CPP4h', 'CPP6h', 'PPO1', 'PPO2', 'I1', 'Iz', 'I2', 'AFp3h', 'AFp4h',...
    'AFF5h', 'AFF6h', 'FFT7h', 'FFC1h', 'FFC2h', 'FFT8h', 'FTT9h', 'FTT7h', 'FCC1h', 'FCC2h', 'FTT8h',...
    'FTT10h', 'TTP7h', 'CCP1h', 'CCP2h', 'TTP8h', 'TPP7h', 'CPP1h', 'CPP2h', 'TPP8h', 'PPO9h', 'PPO5h',...
    'PPO6h', 'PPO10h', 'POO9h', 'POO3h', 'POO4h', 'POO10h', 'OI1h', 'OI2h'};
config.complete_channel_labels = complete_channel_labels;

% Read the performance
performance_human = read_performance_dataset_human(config);
performance_sEEGnal = read_performance_dataset_sEEGnal(config);

% Plot
plot_IC_distribution(config,performance_human,performance_sEEGnal);

% plot_artifacts_distribution(config,performance_human,performance_sEEGnal);

% plot_badchannels_scatter(config,performance_sEEGnal,performance_human)

%%%% FUNCTIONS
% Read all the performance files
function performance = read_performance_dataset_sEEGnal(config)

% Load the dataset
dataset_path = fullfile(config.path.derivatives,'sEEgnal','sEEGnal_dataset.mat');
dummy = load(dataset_path);

% Output structure
modules = {'badchannel_detection', 'artifact_detection'  };
performance = [];
for imodule = 1 : numel(modules)
    performance.times_seconds .(modules{imodule})= nan(numel(dummy.dataset),1);
    performance.memory_bytes.(modules{imodule}) = nan(numel(dummy.dataset),1);
end
performance.badchannels = nan(numel(dummy.dataset),numel(config.complete_channel_labels));
performance.artifacts = cell(numel(dummy.dataset),1);
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
    performance.badchannels(icurrent,:) = ismember(config.complete_channel_labels,current_performance.badchannels.label);
    performance.artifacts{icurrent} = current_performance.artifacts.label;
    
    % Count the number of rejected ICs and store it
    current_ICs_rejected = current_performance.ICs_rejected.IC;
    current_ICs_rejected = current_ICs_rejected(2:end-1);
    ICs_rejected_index = strcmp(current_ICs_rejected,'n/a');
    current_ICs_rejected = current_ICs_rejected(~ICs_rejected_index);
    dummy_count = strjoin(current_ICs_rejected,',');
    dummy_count = strsplit(dummy_count,',');
    performance.ICs_rejected(icurrent) = numel(dummy_count);
    
end

end


function performance = read_performance_dataset_human(config)

% Get the different testers (except sEEGnal)
testers = dir(sprintf('%s/*',config.path.derivatives));
testers = testers(3:end-1);
testers = {testers.name};

% Output structure
modules = {'badchannel_detection', 'artifact_detection'  };
performance = [];
for imodule = 1 : numel(modules)
    performance.times_seconds .(modules{imodule})= nan(20,numel(testers));
    performance.memory_bytes.(modules{imodule}) = nan(20,numel(testers));
end
performance.badchannels = nan(20,numel(config.complete_channel_labels),numel(testers));
performance.artifacts = cell(20,numel(testers));
performance.ICs_rejected = nan(20,numel(testers));

for itester = 1 : numel(testers)
    
    % Load the dataset
    current_tester = testers{itester};
    dataset_path = fullfile(config.path.derivatives,current_tester,...
        sprintf('%s_dataset.mat',current_tester));
    dummy = load(dataset_path);
    
    for icurrent = 1 : numel(dummy.dataset)
        
        % Load performance
        current_dataset = dummy.dataset(icurrent);
        performance_file = fullfile(current_dataset.performance.path,...
            current_dataset.performance.file);
        current_performance = load(performance_file);
        
        % Store times and bytes
        for imodule = 1 : numel(modules)
            performance.times_seconds.(modules{imodule})(icurrent,itester) = current_performance.times_seconds.(modules{imodule});
            performance.memory_bytes.(modules{imodule})(icurrent,itester) = current_performance.memory_bytes.(modules{imodule}) * 0.000001;
        end
        
        % Store number of badchannels and artifacts
        performance.badchannels(icurrent,:,itester) = ismember(config.complete_channel_labels,current_performance.badchannels.label);
        performance.artifacts{icurrent,itester} = current_performance.artifacts.label;
        
        % Count the number of rejected ICs and store it
        performance.ICs_rejected(icurrent,itester) = numel(current_performance.ICs_rejected.IC);
        
    end
    
end

% Average
for imodule = 1 : numel(modules)
    performance.times_seconds.(modules{imodule}) = nanmean(performance.times_seconds.(modules{imodule}),2);
    performance.memory_bytes.(modules{imodule}) = nanmean(performance.memory_bytes.(modules{imodule}),2);
end
performance.ICs_rejected = nanmean(performance.ICs_rejected,2);

end



function plot_IC_distribution(config,performance_human,performance_sEEGnal)

% IC labels
IC_labels = {'clean','noisy'};

% Count each type and estimate the percentage
n_rejected = round(nanmean(performance_human.ICs_rejected));
n_clean = 126 - n_rejected;
count_human = cat(2,ones(1,n_clean),2*ones(1,n_rejected));

n_rejected = round(nanmean(performance_sEEGnal.ICs_rejected));
n_clean = 126 - n_rejected;
count_sEEGnal = cat(2,ones(1,n_clean),2*ones(1,n_rejected));

count = cat(1,count_human,count_sEEGnal);


% Scatter the head with electrodes
fig = figure('WindowState', 'maximized');
imagesc(count)
axis off
title('IC distribution')

% Save the figure
outfile = sprintf('%s/IC_distribution.svg',config.path.figures);
saveas(fig,outfile);
close(fig);


end



function plot_artifacts_distribution(config,performance_human,performance_sEEGnal)

% sEEGnal
% IC labels
artifacst_labels = {'bad_EOG','bad_muscle','bad_jump','bad_other'};

% Count each type and estimate the percentage
count = [];
for itype = 1 : numel(artifacst_labels)
    
    for i = 1 : 20
        
        dummy = sum(strcmp(artifacst_labels{itype},performance_sEEGnal.artifacts{i})); 
        dummy = itype*ones(1,dummy);
        count = cat(2,count,dummy);
        
    end

end

% Scatter the head with electrodes
fig = figure('WindowState', 'maximized');
imagesc(count)
axis off
title('Artifacts distribution sEEGnal')

% Save the figure
outfile = sprintf('%s/artifact_distribution_sEEGnal.svg',config.path.figures);
saveas(fig,outfile);
close(fig);

% Human
% IC labels
artifacst_labels = {'eog','muscle','jump','visual'};

% Count each type and estimate the percentage
count = [];
for itype = 1 : numel(artifacst_labels)
    
    for i = 1 : 20
        
        dummy = sum(strcmp(artifacst_labels{itype},performance_human.artifacts{i})); 
        dummy = itype*ones(1,dummy);
        count = cat(2,count,dummy);
        
    end

end

% Scatter the head with electrodes
fig = figure('WindowState', 'maximized');
imagesc(count)
axis off
title('Artifacts distribution human')

% Save the figure
outfile = sprintf('%s/artifact_distribution_human.svg',config.path.figures);
saveas(fig,outfile);
close(fig);


end



function plot_badchannels_scatter(config,performance_sEEGnal,performance_human)

% Sort the channels based on the number of times is marked as badchannel
percentage = sum(performance_sEEGnal.badchannels,1)/size(performance_sEEGnal.badchannels,1)*100;
[percentage_sorted ,sorted_index] = sort(percentage,'descend');

% Create the figure
fig = figure('WindowState', 'maximized');
hold on

% Plot the channels name in axis X
last = find(percentage_sorted > 0, 1,'last');
xticks(1:last);
dummy_channels = config.complete_channel_labels(sorted_index);
dummy_channels = dummy_channels(1:last);
xticklabels(dummy_channels)
xtickangle(60)

% Add separation lines for plot beauty
line([0 last],[20 20],'Color','black','LineStyle','--')
line([0 last],[10 10],'Color','black','LineStyle','--')

% Divide the channels in four groups and assign a color to each group
thresholds = flipud([0:10:45;10:10:50]');
colors = [0 0 0; 1 0 0; 0.6350,   0.0780, 0.1840;0.8500,0.3250, 0.0980;...
    0.9290,   0.6940, 0.1250; 0.3451    0.8471    0.4392];

% Save the colors to plot later in the head
head_color = nan(numel(config.complete_channel_labels),3);
counter = 1;
for igroup = 1 : size(thresholds,1)
    
    % Get the channels and plot them
    current_threshold = thresholds(igroup,:);
    current_channels_mask = percentage_sorted > current_threshold(1) & percentage_sorted <= current_threshold(2);
    current_x = find(current_channels_mask);
    current_y = percentage_sorted(current_channels_mask); 
    scatter(current_x,current_y,100,colors(igroup,:),'filled');
    
    % Add a line for plot beauty
    for j = 1 : numel(current_x)
        
        line([current_x(j) current_x(j)],[0 current_y(j)],...
            'Color', colors(igroup,:))
        
        % Save the color for later
        head_color(sorted_index(counter),:) = colors(igroup,:);
        counter = counter + 1;
        
    end
    
end

% Assign the color to the "clean" channels
for j = counter : numel(config.complete_channel_labels)
    head_color(sorted_index(j),:) = colors(end,:);
end

% Plot head with colors
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
title('Badchannels sEEGnal')

% Save the figure
outfile = sprintf('%s/badchannels_scatter_sEEGnal.svg',config.path.figures);
saveas(fig,outfile);
close(fig);

% Scatter the head with electrodes
fig = figure('WindowState', 'maximized');
[pos_elec, size_elec] = draw_head(config);
size_elec(:) = 200;
scatter(pos_elec(:,1),pos_elec(:,2),size_elec,head_color,'filled',...
    'MarkerEdgeColor','k')
axis equal off
title('sEEGnal')

% Save the figure
outfile = sprintf('%s/badchannels_head_sEEGnal.svg',config.path.figures);
saveas(fig,outfile);
close(fig);


% HUMANS
% Sort the channels based on the number of times is marked as badchannel
performance_human.badchannels = sum(performance_human.badchannels,3) > 2;
percentage = sum(performance_human.badchannels,1)/size(performance_human.badchannels,1)*100;
[percentage_sorted ,sorted_index] = sort(percentage,'descend');

% Create the figure
fig = figure('WindowState', 'maximized');
hold on

% Plot the channels name in axis X
last = find(percentage_sorted > 0, 1,'last');
xticks(1:last);
dummy_channels = config.complete_channel_labels(sorted_index);
dummy_channels = dummy_channels(1:last);
xticklabels(dummy_channels)
xtickangle(60)

% Add separation lines for plot beauty
line([0 last],[20 20],'Color','black','LineStyle','--')
line([0 last],[10 10],'Color','black','LineStyle','--')

% Divide the channels in four groups and assign a color to each group
thresholds = flipud([0:10:25;10:10:30]');
colors = [0.6350,   0.0780, 0.1840;0.8500,0.3250, 0.0980;...
    0.9290,   0.6940, 0.1250; 0.3451    0.8471    0.4392];

% Save the colors to plot later in the head
head_color = nan(numel(config.complete_channel_labels),3);
counter = 1;
for igroup = 1 : size(thresholds,1)
    
    % Get the channels and plot them
    current_threshold = thresholds(igroup,:);
    current_channels_mask = percentage_sorted > current_threshold(1) & percentage_sorted <= current_threshold(2);
    current_x = find(current_channels_mask);
    current_y = percentage_sorted(current_channels_mask); 
    scatter(current_x,current_y,100,colors(igroup,:),'filled');
    
    % Add a line for plot beauty
    for j = 1 : numel(current_x)
        
        line([current_x(j) current_x(j)],[0 current_y(j)],...
            'Color', colors(igroup,:))
        
        % Save the color for later
        head_color(sorted_index(counter),:) = colors(igroup,:);
        counter = counter + 1;
        
    end
    
end

% Assign the color to the "clean" channels
for j = counter : numel(config.complete_channel_labels)
    head_color(sorted_index(j),:) = colors(end,:);
end

% Plot head with colors
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
title('Badchannels humans')

% Save the figure
outfile = sprintf('%s/badchannels_scatter_humans.svg',config.path.figures);
saveas(fig,outfile);
close(fig);

% Scatter the head with electrodes
fig = figure('WindowState', 'maximized');
[pos_elec, size_elec] = draw_head(config);

size_elec(:) = 200;
scatter(pos_elec(:,1),pos_elec(:,2),size_elec,head_color,'filled',...
    'MarkerEdgeColor','k')
axis equal off
title('Humans')

% Save the figure
outfile = sprintf('%s/badchannels_head_humans.svg',config.path.figures);
saveas(fig,outfile);
close(fig);


end


