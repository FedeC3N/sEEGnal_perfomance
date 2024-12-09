clear
clc
restoredefaultpath

% Add paths
addpath('../shared/')

% Paths
config.path.clean_data = '../../../../databases/LEMON_database/derivatives';
config.path.results = '../../../../results/plv';
if ~exist(config.path.results), mkdir(config.path.results),end

% Get the different testers
testers = {'lemon','sEEGnal'};

% Read the power spectrum
[plv_lemon_dataset,~,channels_lemon_included] = read_plv_dataset(config,'lemon');
[plv_sEEGnal_dataset,f,channels_sEEGnal_included] = read_plv_dataset(config,'sEEGnal');

% Areas
complete_channel_labels = {'Fp1', 'Fp2', 'F7', 'F3', 'Fz', 'F4', 'F8', 'FC5', 'FC1', 'FC2', 'FC6', 'T7', 'C3', 'Cz', 'C4', 'T8', 'CP5', 'CP1', 'CP2', 'CP6', 'AFz', 'P7', 'P3', 'Pz', 'P4', 'P8', 'PO9', 'O1', 'Oz', 'O2', 'PO10', 'AF7', 'AF3', 'AF4', 'AF8', 'F5', 'F1', 'F2', 'F6', 'FT7', 'FC3', 'FC4', 'FT8', 'C5', 'C1', 'C2', 'C6', 'TP7', 'CP3', 'CPz', 'CP4', 'TP8', 'P5', 'P1', 'P2', 'P6', 'PO7', 'PO3', 'POz', 'PO4', 'PO8'};
areas_info = struct('name',{'frontal','temporal_l','temporal_r','parietal_l',...
    'parietal_r','occipital','whole_head'},'channel',[]);
areas_info(1).channel = {'Fp1';'AF7';'AF3'; 'Fpz';'Fp2';'AF8';'AF4';'AFz'};
areas_info(2).channel = {'FT7';'FT9';'T7';'TP7';'TP9'};
areas_info(3).channel = {'FT8';'FT10';'T8';'TP8';'TP10'};
areas_info(4).channel = {'CP5';'CP3';'CP1';'P1';'P3';'P5';'P7';'P9'};
areas_info(5).channel = {'CP6';'CP4';'CP2';'P2';'P4';'P6';'P8';'P10'};
areas_info(6).channel = {'O1';'Iz';'Oz';'O9';'O2';'O10'};
areas_info(7).channel = complete_channel_labels;

% Define the frequency bands
bands_info = struct('name',{'delta', 'theta','alpha','beta','gamma'},...
    'f_limits',{[ 2 4] , [4 8] , [8 12] , [12 20] , [20 45] });

% Stats
% Output struct
results = [];

for iband = 1 : numel(bands_info)
    
    % For saving purposes
    current_band = bands_info(iband).name;
    
    for iarea = 1 : numel(areas_info)
        
        % Get the current band matrix
        current_plv_global_lemon = squeeze(plv_lemon_dataset(:,:,iband,:));
        current_plv_global_sEEGnal = squeeze(plv_sEEGnal_dataset(:,:,iband,:));
        
        % Channels of interest
        desired_channels = areas_info(iarea).channel;
        desired_channels_index = ismember(complete_channel_labels, desired_channels);
        
        % Get the current area matrix
        % The channels of interest vs the rest of channels (just interested
        % in interconnectivity, not intra-connectivity)
        if ~strcmp(areas_info(iarea).name,'whole_head')
            current_plv_global_lemon = current_plv_global_lemon(desired_channels_index,~desired_channels_index,:);
            current_plv_global_sEEGnal = current_plv_global_sEEGnal(desired_channels_index,~desired_channels_index,:);
        end
        
        % Mean connectivity
        current_plv_global_lemon = squeeze(nanmean(nanmean(current_plv_global_lemon,1),2));
        current_plv_global_sEEGnal = squeeze(nanmean(nanmean(current_plv_global_sEEGnal,1),2));
        
        % There are subjects without valid channels. Remove them
        valid = ~isnan(current_plv_global_lemon) & ~isnan(current_plv_global_sEEGnal);
        current_plv_global_lemon = current_plv_global_lemon(valid);
        current_plv_global_sEEGnal = current_plv_global_sEEGnal(valid);
        
        % ttest for the current band -area
        [~,p,~,current_stats] = ttest(current_plv_global_lemon,current_plv_global_sEEGnal);
        
        % Estimate the Effect Size (Cohen's d)
        diff = abs(current_plv_global_lemon - current_plv_global_sEEGnal);
        d = mean(diff)/std(diff);
        
        % mean and stderror
        mean_lemon = mean(current_plv_global_lemon);
        std_mean_error_lemon = std(current_plv_global_lemon)/numel(current_plv_global_lemon);
        mean_sEEGnal = mean(current_plv_global_sEEGnal);
        std_mean_error_sEEGnal = std(current_plv_global_sEEGnal)/numel(current_plv_global_sEEGnal);
        
        % Save
        results.stats.(current_band).(areas_info(iarea).name).test = 'Paired ttest';
        results.stats.(current_band).(areas_info(iarea).name).p = p;
        results.stats.(current_band).(areas_info(iarea).name).stats = current_stats;
        results.stats.(current_band).(areas_info(iarea).name).effect_size_name = 'Cohen d';
        results.stats.(current_band).(areas_info(iarea).name).effect_size = d;
        results.stats.(current_band).(areas_info(iarea).name).mean_lemon = mean_lemon;
        results.stats.(current_band).(areas_info(iarea).name).std_mean_error_lemon = std_mean_error_lemon;
        results.stats.(current_band).(areas_info(iarea).name).mean_sEEGnal = mean_sEEGnal;
        results.stats.(current_band).(areas_info(iarea).name).std_mean_error_sEEGnal = std_mean_error_sEEGnal;
        
        
    end
    
end

% Complete the information
results.areas_info = areas_info;
results.bands_info = bands_info;

% Save the file
outfile = sprintf('%s/plv_results.mat',config.path.results);
save(outfile,'-struct','results');


% Aux functions
function [plv_dataset,bands_info,channels] = read_plv_dataset(config, dataset_name)

% Load the datset
dataset_path = sprintf('%s/%s/%s_dataset.mat',config.path.clean_data,...
    dataset_name,dataset_name);
dummy = load(dataset_path);

for icurrent = 1 : numel(dummy.dataset)
    
    % Load plv
    plv = load(sprintf('%s/%s',dummy.dataset(icurrent).plv.path,...
        dummy.dataset(icurrent).plv.file));
    
    % Save bands_info
    bands_info = plv.bands_info;
    
    % Create the PLV_all struct and the channels info
    if icurrent == 1
        n_sensors = numel(plv.channels_included_index);
        n_bands = numel(bands_info);
        n_subjects = numel(dummy.dataset);
        plv_dataset = nan(n_sensors,n_sensors,...
            n_bands,n_subjects);
        channels = struct('channels_included',[],...
            'channels_included_index',[]);
    end
    
    % Save each PLV
    for iband = 1 : numel(bands_info)
        
        plv_dataset(:,:,iband,icurrent) = plv.(bands_info(iband).name).plv;
        
    end
    
    % Save the channel information
    channels(icurrent).channels_included = plv.channels_included;
    channels(icurrent).channels_included_index = plv.channels_included_index;
    
end

end


% Global connectivity
function test_global_plv(config,plv_SRM_dataset,plv_sEEGnal_dataset,bands_info)

% Output struct
results = [];

for iband = 1 : numel(bands_info)
    
    % Get the current band matrix
    current_plv_global_SRM = squeeze(plv_SRM_dataset(:,:,iband,:));
    current_plv_global_sEEGnal = squeeze(plv_sEEGnal_dataset(:,:,iband,:));
    
    % Get the global connectivity
    current_plv_global_SRM = squeeze(nanmean(nanmean(current_plv_global_SRM,1),2));
    current_plv_global_sEEGnal = squeeze(nanmean(nanmean(current_plv_global_sEEGnal,1),2));
    
    % ttest for the current band
    [~,p,~,current_stats] = ttest(current_plv_global_SRM,current_plv_global_sEEGnal);
    
    % Estimate the Effect Size (Cohen's d)
    diff = abs(current_plv_global_SRM - current_plv_global_sEEGnal);
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
