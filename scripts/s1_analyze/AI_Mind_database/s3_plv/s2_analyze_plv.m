clear
clc
restoredefaultpath

% Add paths
addpath('../shared/')

% Paths
config.path.dataset = '../../../../metadata/AI_Mind_database/dataset';
config.path.stats = '../../../../data/AI_Mind_database/stats';

% Load the whole dataset
load(sprintf('%s/AI_Mind_dataset.mat',config.path.dataset));

% Load the power
[plv_eeg_expert_dataset, ~, channels_SRM] = read_plv_user_dataset(dataset);
[plv_ETL_dataset, bands_info, channels_ETL] = read_plv_dataset(dataset,'etl');

% Areas
complete_channel_labels = {'Fp1', 'Fpz', 'Fp2', 'F7', 'F3', 'Fz', 'F4', 'F8', 'FC5', 'FC1', 'FC2', 'FC6', 'M1', 'T7', 'C3', 'Cz', 'C4', 'T8', 'M2', 'CP5', 'CP1', 'CP2', 'CP6', 'P7', 'P3', 'Pz', 'P4', 'P8', 'POz', 'O1', 'O2', 'AF7', 'AF3', 'AF4', 'AF8', 'F5', 'F1', 'F2', 'F6', 'FC3', 'FCz', 'FC4', 'C5', 'C1', 'C2', 'C6', 'CP3', 'CP4', 'P5', 'P1', 'P2', 'P6', 'F9', 'PO3', 'PO4', 'F10', 'FT7', 'FT8', 'TP7', 'TP8', 'PO7', 'PO8', 'FT9', 'FT10', 'TPP9h', 'TPP10h', 'PO9', 'PO10', 'P9', 'P10', 'AFF1', 'AFz', 'AFF2', 'FFC5h', 'FFC3h', 'FFC4h', 'FFC6h', 'FCC5h', 'FCC3h', 'FCC4h', 'FCC6h', 'CCP5h', 'CCP3h', 'CCP4h', 'CCP6h', 'CPP5h', 'CPP3h', 'CPP4h', 'CPP6h', 'PPO1', 'PPO2', 'I1', 'Iz', 'I2', 'AFp3h', 'AFp4h', 'AFF5h', 'AFF6h', 'FFT7h', 'FFC1h', 'FFC2h', 'FFT8h', 'FTT9h', 'FTT7h', 'FCC1h', 'FCC2h', 'FTT8h', 'FTT10h', 'TTP7h', 'CCP1h', 'CCP2h', 'TTP8h', 'TPP7h', 'CPP1h', 'CPP2h', 'TPP8h', 'PPO9h', 'PPO5h', 'PPO6h', 'PPO10h', 'POO9h', 'POO3h', 'POO4h', 'POO10h', 'OI1h', 'OI2h'}';
areas_info = struct('name',{'frontal','temporal_l','temporal_r','parietal_l',...
    'parietal_r','occipital','whole_head'},'channel',[]);
areas_info(1).channel = {'Fp1';'AF7';'AF3'; 'Fpz';'Fp2';'AF8';'AF4';'AFz'};
areas_info(2).channel = {'FT7';'FT9';'T7';'TP7';'TP9'};
areas_info(3).channel = {'FT8';'FT10';'T8';'TP8';'TP10'};
areas_info(4).channel = {'CP5';'CP3';'CP1';'P1';'P3';'P5';'P7';'P9'};
areas_info(5).channel = {'CP6';'CP4';'CP2';'P2';'P4';'P6';'P8';'P10'};
areas_info(6).channel = {'O1';'Iz';'Oz';'O9';'O2';'O10'};
areas_info(7).channel = complete_channel_labels;

% Stats
% Output struct
results = [];

for iband = 1 : numel(bands_info)
    
    % For saving purposes
    current_band = bands_info(iband).name;
    
    for iarea = 1 : numel(areas_info)
        
        % Get the current band matrix
        current_plv_global_eeg_expert = squeeze(plv_eeg_expert_dataset(:,:,iband,:));
        current_plv_global_ETL = squeeze(plv_ETL_dataset(:,:,iband,:));
        
        % Channels of interest
        desired_channels = areas_info(iarea).channel;
        desired_channels_index = ismember(complete_channel_labels, desired_channels);
        
        % Get the current area matrix
        % The channels of interest vs the rest of channels (just interested
        % in interconnectivity, not intra-connectivity)
        if ~strcmp(areas_info(iarea).name,'whole_head')
            current_plv_global_eeg_expert = current_plv_global_eeg_expert(desired_channels_index,~desired_channels_index,:);
            current_plv_global_ETL = current_plv_global_ETL(desired_channels_index,~desired_channels_index,:);
        end
        
        % Mean connectivity
        current_plv_global_eeg_expert = squeeze(nanmean(nanmean(current_plv_global_eeg_expert,1),2));
        current_plv_global_ETL = squeeze(nanmean(nanmean(current_plv_global_ETL,1),2));
        
        % There are subjects without valid channels. Remove them
        valid = ~isnan(current_plv_global_eeg_expert) & ~isnan(current_plv_global_ETL);
        current_plv_global_eeg_expert = current_plv_global_eeg_expert(valid);
        current_plv_global_ETL = current_plv_global_ETL(valid);
        
        % ttest for the current band -area
        [~,p,~,current_stats] = ttest(current_plv_global_eeg_expert,current_plv_global_ETL);
        
        % Estimate the Effect Size (Cohen's d)
        diff = abs(current_plv_global_eeg_expert - current_plv_global_ETL);
        d = mean(diff)/std(diff);
        
        % mean and stderror
        mean_eeg_expert = mean(current_plv_global_eeg_expert);
        std_mean_error_eeg_expert = std(current_plv_global_eeg_expert)/numel(current_plv_global_eeg_expert);
        mean_ETL = mean(current_plv_global_ETL);
        std_mean_error_ETL = std(current_plv_global_ETL)/numel(current_plv_global_ETL);
        
        % Save
        results.stats.(current_band).(areas_info(iarea).name).test = 'Paired ttest';
        results.stats.(current_band).(areas_info(iarea).name).p = p;
        results.stats.(current_band).(areas_info(iarea).name).stats = current_stats;
        results.stats.(current_band).(areas_info(iarea).name).effect_size_name = 'Cohen d';
        results.stats.(current_band).(areas_info(iarea).name).effect_size = d;
        results.stats.(current_band).(areas_info(iarea).name).mean_eeg_expert = mean_eeg_expert;
        results.stats.(current_band).(areas_info(iarea).name).std_mean_error_eeg_expert = std_mean_error_eeg_expert;
        results.stats.(current_band).(areas_info(iarea).name).mean_ETL = mean_ETL;
        results.stats.(current_band).(areas_info(iarea).name).std_mean_error_ETL = std_mean_error_ETL;
        
        
    end
    
end

% Complete the information
results.areas_info = areas_info;
results.bands_info = bands_info;

% Save the file
outfile = sprintf('%s/plv_stats.mat',config.path.stats);
save(outfile,'-struct','results');


% Aux functions
function [plv_dataset,bands_info,channels] = read_plv_dataset(dataset, desired_dataset)

% subject of interest
current_dataset_index = ismember({dataset.origin},desired_dataset);
current_dataset = dataset(current_dataset_index);

for idataset = 1 : numel(current_dataset)
    
    % Load pow
    plv = load(sprintf('%s/%s',current_dataset(idataset).plv.path,...
        current_dataset(idataset).plv.file));
    
    % Save bands_info
    bands_info = plv.bands_info;
    
    % Create the PLV_all struct and the channels info
    if idataset == 1
        n_sensors = numel(plv.channels_included_index);
        n_bands = numel(bands_info);
        n_subjects = numel(current_dataset);
        plv_dataset = nan(n_sensors,n_sensors,...
            n_bands,n_subjects);
        channels = struct('channels_included',[],...
            'channels_included_index',[]);
    end
    
    % Save each PLV
    for iband = 1 : numel(bands_info)
        
        plv_dataset(:,:,iband,idataset) = plv.(bands_info(iband).name).plv;
        
    end
    
    % Save the channel information
    channels(idataset).channels_included = plv.channels_included;
    channels(idataset).channels_included_index = plv.channels_included_index;
    
end

end


function [plv_eeg_expert_dataset,bands_info,channels_eeg_expert] = read_plv_user_dataset(dataset)


users = {dataset.origin};
users = unique(users(~ismember(users,'etl')));
for iuser = 1 : numel(users)
    
    [current_plv_eeg_expert_dataset,bands_info,channels_eeg_expert] = ...
        read_plv_dataset(dataset,users{iuser});
    
    if iuser == 1
        n_sensors = size(current_plv_eeg_expert_dataset,1);
        n_bands = size(current_plv_eeg_expert_dataset,3);
        n_subjects = size(current_plv_eeg_expert_dataset,4);
        plv_eeg_expert_dataset = nan(n_sensors,n_sensors,...
            n_bands,n_subjects,numel(users));
    end
    plv_eeg_expert_dataset(:,:,:,:,iuser) = current_plv_eeg_expert_dataset;
    
end

plv_eeg_expert_dataset = nanmean(plv_eeg_expert_dataset,5);

end

% Global connectivity
function test_global_plv(config,plv_SRM_dataset,plv_ETL_dataset,bands_info)

% Output struct
results = [];

for iband = 1 : numel(bands_info)
    
    % Get the current band matrix
    current_plv_global_SRM = squeeze(plv_SRM_dataset(:,:,iband,:));
    current_plv_global_ETL = squeeze(plv_ETL_dataset(:,:,iband,:));
    
    % Get the global connectivity
    current_plv_global_SRM = squeeze(nanmean(nanmean(current_plv_global_SRM,1),2));
    current_plv_global_ETL = squeeze(nanmean(nanmean(current_plv_global_ETL,1),2));
    
    % ttest for the current band
    [~,p,~,current_stats] = ttest(current_plv_global_SRM,current_plv_global_ETL);
    
    % Estimate the Effect Size (Cohen's d)
    diff = abs(current_plv_global_SRM - current_plv_global_ETL);
    d = mean(diff)/std(diff);
    
    % Save
    current_band = bands_info(iband).name;
    results.stats.(current_band).test = 'Paired ttest';
    results.stats.(current_band).p = p;
    results.stats.(current_band).stats = current_stats;
    results.stats.(current_band).effect_size_name = 'Cohen d';
    results.stats.(current_band).effect_size = d;
    
    
end

% Complete the information
results.bands_info = bands_info;

% Save the file
outfile = sprintf('%s/plv_stats.mat',config.path.stats);
save(outfile,'-struct','results');

end
