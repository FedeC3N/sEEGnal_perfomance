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

% To define later the pow matrix
complete_channel_labels = {'Fp1', 'Fpz', 'Fp2', 'F7', 'F3', 'Fz', 'F4', 'F8', 'FC5', 'FC1', 'FC2', 'FC6', 'M1', 'T7', 'C3', 'Cz', 'C4', 'T8', 'M2', 'CP5', 'CP1', 'CP2', 'CP6', 'P7', 'P3', 'Pz', 'P4', 'P8', 'POz', 'O1', 'O2', 'AF7', 'AF3', 'AF4', 'AF8', 'F5', 'F1', 'F2', 'F6', 'FC3', 'FCz', 'FC4', 'C5', 'C1', 'C2', 'C6', 'CP3', 'CP4', 'P5', 'P1', 'P2', 'P6', 'F9', 'PO3', 'PO4', 'F10', 'FT7', 'FT8', 'TP7', 'TP8', 'PO7', 'PO8', 'FT9', 'FT10', 'TPP9h', 'TPP10h', 'PO9', 'PO10', 'P9', 'P10', 'AFF1', 'AFz', 'AFF2', 'FFC5h', 'FFC3h', 'FFC4h', 'FFC6h', 'FCC5h', 'FCC3h', 'FCC4h', 'FCC6h', 'CCP5h', 'CCP3h', 'CCP4h', 'CCP6h', 'CPP5h', 'CPP3h', 'CPP4h', 'CPP6h', 'PPO1', 'PPO2', 'I1', 'Iz', 'I2', 'AFp3h', 'AFp4h', 'AFF5h', 'AFF6h', 'FFT7h', 'FFC1h', 'FFC2h', 'FFT8h', 'FTT9h', 'FTT7h', 'FCC1h', 'FCC2h', 'FTT8h', 'FTT10h', 'TTP7h', 'CCP1h', 'CCP2h', 'TTP8h', 'TPP7h', 'CPP1h', 'CPP2h', 'TPP8h', 'PPO9h', 'PPO5h', 'PPO6h', 'PPO10h', 'POO9h', 'POO3h', 'POO4h', 'POO10h', 'OI1h', 'OI2h'}';

% Load the power
[pow_ETL_dataset_norm,f,channels_ETL] = read_pow_dataset(dataset,'etl');
[pow_eeg_expert_dataset_norm,f,channels_eeg_expert] = read_pow_user_dataset(dataset);


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
        avg_band_pow_eeg_expert = nan(1,size(pow_eeg_expert_dataset_norm,3));
        avg_band_pow_ETL = nan(1,size(pow_ETL_dataset_norm,3));
        
        % Selects the channel for the current subject
        desired_channels = areas_info(iarea).channel;
        desired_channels_index = ismember(complete_channel_labels, desired_channels);
        
        % and estimate the average power
        for isubject = 1 : size(pow_eeg_expert_dataset_norm,3)
            
            % EEG experts dataset
            current_channels = desired_channels_index & channels_eeg_expert(isubject).channels_included_index;
            current_pow_norm = pow_eeg_expert_dataset_norm(current_channels,current_f,isubject);
            current_pow_norm = mean(mean(current_pow_norm,1),2);
            
            avg_band_pow_eeg_expert(isubject) = current_pow_norm;
        end
        
        % ETL dataset
        for isubject = 1 : size(pow_ETL_dataset_norm,3)
            
            current_channels = desired_channels_index & channels_ETL(isubject).channels_included_index;
            current_pow_norm = pow_ETL_dataset_norm(current_channels,current_f,isubject);
            current_pow_norm = mean(mean(current_pow_norm,1),2);
            
            avg_band_pow_ETL(isubject) = current_pow_norm;
            
        end
        
        % There are subjects without valid channels. Remove them
        valid = ~isnan(avg_band_pow_eeg_expert) & ~isnan(avg_band_pow_ETL);
        avg_band_pow_eeg_expert = avg_band_pow_eeg_expert(valid);
        avg_band_pow_ETL = avg_band_pow_ETL(valid);
        
        % ttest for the current_area - band
        [~,p,~,current_stats] = ttest(avg_band_pow_eeg_expert,avg_band_pow_ETL);
        
        % Estimate the Effect Size (Cohen's d)
        diff = abs(avg_band_pow_eeg_expert - avg_band_pow_ETL);
        d = mean(diff)/std(diff);
        
        % mean and stderror
        mean_eeg_expert = mean(avg_band_pow_eeg_expert);
        std_mean_error_eeg_expert = std(avg_band_pow_eeg_expert)/numel(avg_band_pow_eeg_expert);
        mean_ETL = mean(avg_band_pow_ETL);
        std_mean_error_ETL = std(avg_band_pow_ETL)/numel(avg_band_pow_ETL);
        
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
outfile = sprintf('%s/pow_stats.mat',config.path.stats);
save(outfile,'-struct','results');



% Functions
function [pow_dataset_norm,f,channels] = read_pow_dataset(dataset, desired_dataset)

% subject of interest
current_dataset_index = ismember({dataset.origin},desired_dataset);
current_dataset = dataset(current_dataset_index);

for icurrent = 1 : numel(current_dataset)
    
    % Load pow
    pow = load(sprintf('%s/%s',current_dataset(icurrent).pow.path,...
        current_dataset(icurrent).pow.file));
    
    % Save f
    f = pow.f;
    
    % Normalize
    current_pow = mean(pow.pow_spectrum,3);
    scaling_factor = 1./nansum(nansum(current_pow,2));
    scaling_factor = repmat(scaling_factor,size(current_pow));
    current_pow_norm = scaling_factor .* current_pow;
    
    % Add to the all matrix
    if icurrent == 1
        pow_dataset_norm = nan(64,numel(f),numel(current_dataset));
        channels = struct('channels_included',[],...
            'channels_included_index',[]);
    end
    pow_dataset_norm(1:size(current_pow_norm,1),:,icurrent) = current_pow_norm;
    channels(icurrent).channels_included = pow.channels_included;
    channels(icurrent).channels_included_index = pow.channels_included_index;
    
end

end


function [pow_eeg_expert_dataset_norm,f,channels_eeg_expert] = read_pow_user_dataset(dataset)

users = {dataset.origin};
users = unique(users(~ismember(users,'etl')));
for iuser = 1 : numel(users)
    
    [current_pow_eeg_expert_dataset_norm,f,channels_eeg_expert] = ...
        read_pow_dataset(dataset,users{iuser});
    
    if iuser == 1
        pow_eeg_expert_dataset_norm = nan([size(current_pow_eeg_expert_dataset_norm),numel(users)]);
    end
    pow_eeg_expert_dataset_norm(:,:,:,iuser) = current_pow_eeg_expert_dataset_norm;
    
end

pow_eeg_expert_dataset_norm = nanmean(pow_eeg_expert_dataset_norm,4);

end
