clear
clc
restoredefaultpath

% Paths
config.path.dataset = '../../../../data/SRM_database/dataset';
config.path.stats = '../../../../data/SRM_database/stats';

% Create output path
if ~exist(config.path.stats), mkdir(config.path.stats), end

% Load the whole dataset
load(sprintf('%s/SRM_dataset.mat',config.path.dataset));

% Load the IAF
iaf_SRM_dataset = read_iaf_dataset(dataset,'SRM_database');
iaf_ETL_dataset = read_iaf_dataset(dataset,'ETL_database');

% Measuser
measures = {'iaf', 'iaf_amp'};

% Remove nans
index_SRM = isnan([iaf_SRM_dataset.iaf]);
index_ETL = isnan([iaf_ETL_dataset.iaf]);
index_nan = index_SRM | index_ETL;
iaf_SRM_dataset(index_nan) = [];
iaf_ETL_dataset(index_nan) = [];

% Variable for the stats
results = [];

for imeasure = 1 : numel(measures)
    
    current_measure = measures{imeasure};
    
    % ttest for iaf
    [~,p,~,current_stats] = ttest([iaf_SRM_dataset.(current_measure)],[iaf_ETL_dataset.(current_measure)]);
    
    % Estimate the Effect Size (Cohen's d)
    diff = abs([iaf_SRM_dataset.(current_measure)] - [iaf_ETL_dataset.(current_measure)]);
    d = mean(diff)/std(diff);
    
    % mean and stderror
    mean_SRM = mean([iaf_SRM_dataset.(current_measure)]);
    std_mean_error_SRM = std([iaf_SRM_dataset.(current_measure)])/numel([iaf_SRM_dataset.(current_measure)]);
    mean_ETL = mean([iaf_ETL_dataset.(current_measure)]);
    std_mean_error_ETL = std([iaf_ETL_dataset.(current_measure)])/numel([iaf_ETL_dataset.(current_measure)]);
    
    % Save
    results.stats.(current_measure).test = 'Paired ttest';
    results.stats.(current_measure).p = p;
    results.stats.(current_measure).stats = current_stats;
    results.stats.(current_measure).effect_size_name = 'Cohen d';
    results.stats.(current_measure).effect_size = d;
    results.stats.(current_measure).mean_SRM = mean_SRM;
    results.stats.(current_measure).std_mean_error_SRM = std_mean_error_SRM;
    results.stats.(current_measure).mean_ETL = mean_ETL;
    results.stats.(current_measure).std_mean_error_ETL = std_mean_error_ETL;
end

% Save the file
outfile = sprintf('%s/iaf_stats.mat',config.path.stats);
save(outfile,'-struct','results');


% Functions
function iaf_dataset = read_iaf_dataset(dataset, desired_dataset)

% subject of interest
current_dataset_index = ismember({dataset.database},desired_dataset);
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
