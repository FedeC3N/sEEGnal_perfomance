clear
clc
restoredefaultpath

% Paths
config.path.dataset = '../../../../metadata/AI_Mind_database/dataset';
config.path.stats = '../../../../data/AI_Mind_database/stats';

% Create output path
if ~exist(config.path.stats), mkdir(config.path.stats), end

% Load the whole dataset
load(sprintf('%s/AI_Mind_dataset.mat',config.path.dataset));

% Load the IAF
iaf_eeg_expert_dataset = read_iaf_eeg_expert_dataset(dataset);
iaf_ETL_dataset = read_iaf_dataset(dataset,'etl');

% Measures
measures = {'iaf', 'iaf_amp'};

% Remove nans
index_eeg_expert = isnan([iaf_eeg_expert_dataset.iaf]);
index_ETL = isnan([iaf_ETL_dataset.iaf]);
index_nan = index_eeg_expert | index_ETL;
iaf_eeg_expert_dataset(index_nan) = [];
iaf_ETL_dataset(index_nan) = [];

% Variable for the stats
results = [];

for imeasure = 1 : numel(measures)
    
    current_measure = measures{imeasure};
    
    % ttest for iaf
    [~,p,~,current_stats] = ttest([iaf_eeg_expert_dataset.(current_measure)],[iaf_ETL_dataset.(current_measure)]);
    
    % Estimate the Effect Size (Cohen's d)
    diff = abs([iaf_eeg_expert_dataset.(current_measure)] - [iaf_ETL_dataset.(current_measure)]);
    d = mean(diff)/std(diff);
    
    % mean and stderror
    mean_eeg_expert = mean([iaf_eeg_expert_dataset.(current_measure)]);
    std_mean_error_eeg_expert = std([iaf_eeg_expert_dataset.(current_measure)])/numel([iaf_eeg_expert_dataset.(current_measure)]);
    mean_ETL = mean([iaf_ETL_dataset.(current_measure)]);
    std_mean_error_ETL = std([iaf_ETL_dataset.(current_measure)])/numel([iaf_ETL_dataset.(current_measure)]);
    
    % Save
    results.stats.(current_measure).test = 'Paired ttest';
    results.stats.(current_measure).p = p;
    results.stats.(current_measure).stats = current_stats;
    results.stats.(current_measure).effect_size_name = 'Cohen d';
    results.stats.(current_measure).effect_size = d;
    results.stats.(current_measure).mean_eeg_expert = mean_eeg_expert;
    results.stats.(current_measure).std_mean_error_eeg_expert = std_mean_error_eeg_expert;
    results.stats.(current_measure).mean_ETL = mean_ETL;
    results.stats.(current_measure).std_mean_error_ETL = std_mean_error_ETL;
end

% Save the file
outfile = sprintf('%s/iaf_stats.mat',config.path.stats);
save(outfile,'-struct','results');


% Functions
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


function iaf_eeg_expert_dataset = read_iaf_eeg_expert_dataset(dataset)

users = {dataset.origin};
users = unique(users(~ismember(users,'etl')));
for iuser = 1 : numel(users)
    
    current_iaf_eeg_expert_dataset = ...
        read_iaf_dataset(dataset,users{iuser});
    
    % Little bit complicated struct
    % For the first eeg_expert, store the result
    if iuser == 1
        iaf_eeg_expert_dataset = struct('iaf',[],'iaf_amp',[]);
        
        for isubject = 1 : numel(current_iaf_eeg_expert_dataset)
            iaf_eeg_expert_dataset(isubject).iaf = current_iaf_eeg_expert_dataset(isubject).iaf;
            iaf_eeg_expert_dataset(isubject).iaf_amp = current_iaf_eeg_expert_dataset(isubject).iaf_amp;
        end
    % For the rest, sum them
    else
        for isubject = 1 : numel(current_iaf_eeg_expert_dataset)
            iaf_eeg_expert_dataset(isubject).iaf = iaf_eeg_expert_dataset(isubject).iaf + current_iaf_eeg_expert_dataset(isubject).iaf;
            iaf_eeg_expert_dataset(isubject).iaf_amp = iaf_eeg_expert_dataset(isubject).iaf_amp + current_iaf_eeg_expert_dataset(isubject).iaf_amp;
        end
    end
    
end

% Finally, obtain tha average dividing by the number of eeg_experts
for isubject = 1 : numel(current_iaf_eeg_expert_dataset)
    iaf_eeg_expert_dataset(isubject).iaf = iaf_eeg_expert_dataset(isubject).iaf/numel(users);
    iaf_eeg_expert_dataset(isubject).iaf_amp = iaf_eeg_expert_dataset(isubject).iaf_amp/numel(users);
end

end
