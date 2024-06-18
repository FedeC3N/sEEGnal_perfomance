clear
clc
restoredefaultpath

% Add paths
addpath('../shared/')

% Paths
config.path.dataset = '../../../../data/SRM_database/dataset';
config.path.stats = '../../../../data/SRM_database/stats';

% Load the whole dataset
load(sprintf('%s/SRM_dataset.mat',config.path.dataset));

% Load the power
[plv_SRM_dataset, ~, channels_SRM] = read_plv_dataset(dataset,'SRM_database');
[plv_ETL_dataset, bands_info, channels_ETL] = read_plv_dataset(dataset,'ETL_database');

% Areas
complete_channel_labels = {'Fp1';'AF7';'AF3';'F1';'F3';'F5';'F7';'FT7';'FC5';'FC3';'FC1';'C1';'C3';'C5';'T7';'TP7';'CP5';'CP3';'CP1';'P1';'P3';'P5';'P7';'P9';'PO7';'PO3';'O1';'Iz';'Oz';'POz';'Pz';'CPz';'Fpz';'Fp2';'AF8';'AF4';'AFz';'Fz';'F2';'F4';'F6';'F8';'FT8';'FC6';'FC4';'FC2';'FCz';'Cz';'C2';'C4';'C6';'T8';'TP8';'CP6';'CP4';'CP2';'P2';'P4';'P6';'P8';'P10';'PO8';'PO4';'O2'};
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
        current_plv_global_SRM = squeeze(plv_SRM_dataset(:,:,iband,:));
        current_plv_global_ETL = squeeze(plv_ETL_dataset(:,:,iband,:));
        
        % Channels of interest
        desired_channels = areas_info(iarea).channel;
        desired_channels_index = ismember(complete_channel_labels, desired_channels);

        % Get the current area matrix
        % The channels of interest vs the rest of channels (just interested
        % in interconnectivity, not intra-connectivity)
        if ~strcmp(areas_info(iarea).name,'whole_head')
            current_plv_global_SRM = current_plv_global_SRM(desired_channels_index,~desired_channels_index,:);
            current_plv_global_ETL = current_plv_global_ETL(desired_channels_index,~desired_channels_index,:);
        end
        
        % Mean connectivity
        current_plv_global_SRM = squeeze(nanmean(nanmean(current_plv_global_SRM,1),2));
        current_plv_global_ETL = squeeze(nanmean(nanmean(current_plv_global_ETL,1),2));
    
        % There are subjects without valid channels. Remove them
        valid = ~isnan(current_plv_global_SRM) & ~isnan(current_plv_global_ETL);
        current_plv_global_SRM = current_plv_global_SRM(valid);
        current_plv_global_ETL = current_plv_global_ETL(valid);
        
        % ttest for the current band -area
        [~,p,~,current_stats] = ttest(current_plv_global_SRM,current_plv_global_ETL);
    
        % Estimate the Effect Size (Cohen's d)
        diff = abs(current_plv_global_SRM - current_plv_global_ETL);
        d = mean(diff)/std(diff);
        
        % mean and stderror
        mean_SRM = mean(current_plv_global_SRM);
        std_mean_error_SRM = std(current_plv_global_SRM)/numel(current_plv_global_SRM);
        mean_ETL = mean(current_plv_global_ETL);
        std_mean_error_ETL = std(current_plv_global_ETL)/numel(current_plv_global_ETL);
        
        % Save
        results.stats.(current_band).(areas_info(iarea).name).test = 'Paired ttest';
        results.stats.(current_band).(areas_info(iarea).name).p = p;
        results.stats.(current_band).(areas_info(iarea).name).stats = current_stats;
        results.stats.(current_band).(areas_info(iarea).name).effect_size_name = 'Cohen d';
        results.stats.(current_band).(areas_info(iarea).name).effect_size = d;
        results.stats.(current_band).(areas_info(iarea).name).mean_SRM = mean_SRM;
        results.stats.(current_band).(areas_info(iarea).name).std_mean_error_SRM = std_mean_error_SRM;
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
current_dataset_index = ismember({dataset.database},desired_dataset);
current_dataset = dataset(current_dataset_index);



for idataset = 1 : numel(current_dataset)
    
    % Load pow
    plv = load(sprintf('%s/%s',current_dataset(idataset).plv.path,...
        current_dataset(idataset).plv.file));
    
    % Save bands_info
    bands_info = plv.bands_info;
    
    % Create the PLV_all struct and the channels info
    if idataset == 1
        plv_dataset = nan(64,64,numel(bands_info),numel(current_dataset));
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
