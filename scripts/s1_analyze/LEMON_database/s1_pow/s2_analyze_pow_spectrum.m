clear
clc
restoredefaultpath

% Paths
config.path.clean_data = '../../../../databases/LEMON_database/derivatives';
config.path.results = '../../../../results/pow';

if ~exist(config.path.results),mkdir(config.path.results),end

% To define later the pow matrix
complete_channel_labels = {'Fp1', 'Fp2', 'F7', 'F3', 'Fz', 'F4', 'F8', 'FC5', 'FC1', 'FC2', 'FC6', 'T7', 'C3', 'Cz', 'C4', 'T8', 'CP5', 'CP1', 'CP2', 'CP6', 'AFz', 'P7', 'P3', 'Pz', 'P4', 'P8', 'PO9', 'O1', 'Oz', 'O2', 'PO10', 'AF7', 'AF3', 'AF4', 'AF8', 'F5', 'F1', 'F2', 'F6', 'FT7', 'FC3', 'FC4', 'FT8', 'C5', 'C1', 'C2', 'C6', 'TP7', 'CP3', 'CPz', 'CP4', 'TP8', 'P5', 'P1', 'P2', 'P6', 'PO7', 'PO3', 'POz', 'PO4', 'PO8'};

% Get the different testers
testers = {'lemon','sEEGnal'};

% Read the power spectrum
[pow_lemon_dataset_norm,~,channels_lemon_included] = read_pow_dataset(config,'lemon');
[pow_sEEGnal_dataset_norm,f,channels_sEEGnal_included] = read_pow_dataset(config,'sEEGnal');

% Per band and area
% Define the frequency bands
bands_info = struct('name',{'delta', 'theta','alpha','beta','gamma', 'broadband'},...
    'f_limits',{[ 2 4] , [4 8] , [8 12] , [12 20] , [20 45] , [2 45]},...
    'f_limits_index',[],'f_original',f');

% Areas
areas_info = struct('name',{'frontal','temporal_l','temporal_r','parietal_l',...
    'parietal_r','occipital','whole_head'},'channel',[]);
areas_info(1).channel = {'Fp1';'AF7';'AF3'; 'Fpz';'Fp2';'AF8';'AF4';'AFz'};
areas_info(2).channel = {'FT7';'FT9';'T7';'TP7';'TP9'};
areas_info(3).channel = {'FT8';'FT10';'T8';'TP8';'TP10'};
areas_info(4).channel = {'CP5';'CP3';'CP1';'P1';'P3';'P5';'P7';'P9'};
areas_info(5).channel = {'CP6';'CP4';'CP2';'P2';'P4';'P6';'P8';'P10'};
areas_info(6).channel = {'O1';'Iz';'Oz';'O9';'O2';'O10'};
areas_info(7).channel = complete_channel_labels;

% Variable for the stats
results = [];

for iband = 1 : numel(bands_info)
    
    % Obtain the current freqencies of interest
    current_band = bands_info(iband).name;
    current_band_f_limits = bands_info(iband).f_limits;
    current_f = f>current_band_f_limits(1) & f<current_band_f_limits(2);
    current_f = current_f';
    bands_info(iband).f_limits_index = current_f;
    
    % For each area
    for iarea = 1 : numel (areas_info)
        
        % Save space for the comparison
        avg_band_pow_lemon = nan(1,size(pow_lemon_dataset_norm,3));
        avg_band_pow_sEEGnal = nan(1,size(pow_sEEGnal_dataset_norm,3));
        
        % Selects the channel for the current subject
        desired_channels = areas_info(iarea).channel;
        desired_channels_index = ismember(complete_channel_labels, desired_channels);
        
        % and estimate the average power
        for isubject = 1 : size(pow_lemon_dataset_norm,3)
            
            % EEG experts dataset
            current_channels = desired_channels_index & channels_lemon_included(isubject).channels_included_index;
            current_pow_norm = pow_lemon_dataset_norm(current_channels,current_f,isubject);
            current_pow_norm = mean(mean(current_pow_norm,1),2);
            
            avg_band_pow_lemon(isubject) = current_pow_norm;
        end
        
        % ETL dataset
        for isubject = 1 : size(pow_sEEGnal_dataset_norm,3)
            
            current_channels = desired_channels_index & channels_sEEGnal_included(isubject).channels_included_index;
            current_pow_norm = pow_sEEGnal_dataset_norm(current_channels,current_f,isubject);
            current_pow_norm = mean(mean(current_pow_norm,1),2);
            
            avg_band_pow_sEEGnal(isubject) = current_pow_norm;
            
        end
        
        % There are subjects without valid channels. Remove them
        valid = ~isnan(avg_band_pow_lemon) & ~isnan(avg_band_pow_sEEGnal);
        avg_band_pow_lemon = avg_band_pow_lemon(valid);
        avg_band_pow_sEEGnal = avg_band_pow_sEEGnal(valid);
        
        % ttest for the current_area - band
        [~,p,~,current_stats] = ttest(avg_band_pow_lemon,avg_band_pow_sEEGnal);
        
        % Estimate the Effect Size (Cohen's d)
        diff = abs(avg_band_pow_lemon - avg_band_pow_sEEGnal);
        d = mean(diff)/std(diff);
        
        % mean and stderror
        mean_lemon = mean(avg_band_pow_lemon);
        std_mean_error_lemon = std(avg_band_pow_lemon)/numel(avg_band_pow_lemon);
        mean_sEEGnal = mean(avg_band_pow_sEEGnal);
        std_mean_error_sEEGnal = std(avg_band_pow_sEEGnal)/numel(avg_band_pow_sEEGnal);
        
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
outfile = sprintf('%s/pow_results.mat',config.path.results);
save(outfile,'-struct','results');


% Functions
function [pow_dataset_norm,f,channels] = read_pow_dataset(config,dataset_name)

% Load the datset
dataset_path = sprintf('%s/%s/%s_dataset.mat',config.path.clean_data,...
    dataset_name,dataset_name);
dummy = load(dataset_path);

for icurrent = 1 : numel(dummy.dataset)
    
    % Load pow
    current_dataset = dummy.dataset(icurrent);
    pow_file = sprintf('%s/%s',current_dataset.pow.path,...
        current_dataset.pow.file);
    pow = load(pow_file);
    
    % Save f
    f = pow.f;
    
    % Normalize
    current_pow = mean(pow.pow_spectrum,3);
    scaling_factor = 1./nansum(nansum(current_pow,2));
    scaling_factor = repmat(scaling_factor,size(current_pow));
    current_pow_norm = scaling_factor .* current_pow;
    
    % Add to the all matrix
    if icurrent == 1
        pow_dataset_norm = nan(numel(pow.complete_channel_labels),numel(f),numel(dataset));
        channels = struct('channels_included',[],...
            'channels_included_index',[]);
    end
    pow_dataset_norm(1:size(current_pow_norm,1),:,icurrent) = current_pow_norm;
    channels(icurrent).channels_included = pow.channels_included;
    channels(icurrent).channels_included_index = pow.channels_included_index;
    
end

end

