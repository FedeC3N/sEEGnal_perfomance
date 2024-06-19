clear
clc
close all
restoredefaultpath

% Add paths
addpath('../shared/')

% Paths
config.path.dataset = '../../../../data/SRM_database/dataset';
config.path.stats = '../../../../data/SRM_database/stats';

% Load the whole dataset
load(sprintf('%s/SRM_dataset.mat',config.path.dataset));

% Load the IAF
iaf_SRM_dataset = read_iaf_dataset(dataset,'SRM_database');
iaf_ETL_dataset = read_iaf_dataset(dataset,'ETL_database');

% Measures
measures = {'iaf', 'iaf_amp'};

% Plot difference in channels
plot_diff_channels(iaf_SRM_dataset,iaf_ETL_dataset,measures)

% Aux functions
function iaf_dataset = read_iaf_dataset(dataset, desired_dataset)

% subject of interest
current_dataset_index = ismember({dataset.origin},desired_dataset);
current_dataset = dataset(current_dataset_index);

for icurrent = 1 : numel(current_dataset)
    
    % Load pow
    iaf = load(sprintf('%s/%s',current_dataset(icurrent).iaf.path,...
        current_dataset(icurrent).iaf.file));
    
    % Add to the all matrix
    if icurrent == 1
        iaf_dataset = struct('iaf',[],'iaf_amp',[]);
    end
    iaf_dataset(icurrent).iaf = iaf.iaf;
    iaf_dataset(icurrent).iaf_amp =iaf.iaf_amp;
    
end

end


function plot_diff_channels(iaf_SRM_dataset,iaf_ETL_dataset,measures)

figure('WindowState', 'maximized');
hold on

% For each band
for imeasure = 1 : numel(measures)
    
    % Select the current measure
    current_iaf_eeg_expert = [iaf_SRM_dataset(:).(measures{imeasure})];
    current_iaf_ETL = [iaf_ETL_dataset(:).(measures{imeasure})];
    
    % Get the difference
    iaf_diff_vector = current_iaf_eeg_expert - current_iaf_ETL;
    x_vector = imeasure * ones(numel(iaf_diff_vector),1);
    
    
    % Plot
    sw = swarmchart(x_vector,iaf_diff_vector,'filled','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5);
    bx = boxchart(x_vector,iaf_diff_vector,'BoxFaceColor',sw.CData ./ 1.2,'WhiskerLineColor',sw.CData ./ 1.2,...
        'MarkerStyle','none','BoxWidth',sw.XJitterWidth);
    
    
end

xlim([0 numel(measures) + 1])
xticks(1:numel(measures))
xticklabels({'IAF', 'IAF Amplitude'})
title('Average IAF difference')
set(gca,'TickLabelInterpreter','none')



end






