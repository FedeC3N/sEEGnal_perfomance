clear
clc
close all
restoredefaultpath

% Add paths
addpath('../shared/')

% Paths
config.path.dataset = '../../../../metadata/AI_Mind_database/dataset';
config.path.stats = '../../../../results';

% Load the whole dataset
load(sprintf('%s/AI_Mind_dataset.mat',config.path.dataset));

% Load the power
[plv_eeg_expert_dataset, ~, channels_eeg_expert] = read_plv_user_dataset(dataset);
[plv_ETL_dataset, bands_info, channels_ETL] = read_plv_dataset(dataset,'etl');

% Plot difference in channels
plot_diff_channels(plv_eeg_expert_dataset,plv_ETL_dataset,bands_info)

% Plot correlation of each pair of channels
plot_corr_channels(plv_eeg_expert_dataset,plv_ETL_dataset,bands_info, channels_ETL)

% Plot correlation of PLV matrix for each subject
plot_corr_plv_subjects(plv_eeg_expert_dataset,plv_ETL_dataset,bands_info, channels_eeg_expert,channels_ETL)

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



function plot_diff_channels(plv_eeg_expert_dataset,plv_ETL_dataset,bands_info)

figure('WindowState', 'maximized');
hold on

% For each band
for iband = 1 : numel(bands_info)
    
    % Select the band data
    current_plv_global_eeg_expert = squeeze(plv_eeg_expert_dataset(:,:,iband,:));
    current_plv_global_ETL = squeeze(plv_ETL_dataset(:,:,iband,:));
    
    % Get the difference of each channel and average all the differences
    plv_diff = current_plv_global_eeg_expert - current_plv_global_ETL;
    plv_diff_avg = nanmean(plv_diff,3);
    
    % Symmetrical. Get only upper matrix
    plv_diff_vector = triu(ones(size(plv_diff_avg)),1) > 0 ;
    plv_diff_vector = plv_diff_avg(plv_diff_vector);
    x_vector = iband * ones(numel(plv_diff_vector),1);
    
    % Express the difference as percentage
    plv_diff_vector_perc = (abs(plv_diff_vector)/2)*100;
    
    % Plot
    sw = swarmchart(x_vector,plv_diff_vector_perc,'filled','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5);
    bx = boxchart(x_vector,plv_diff_vector_perc,'BoxFaceColor',sw.CData ./ 1.2,'WhiskerLineColor',sw.CData ./ 1.2,...
        'MarkerStyle','none','BoxWidth',sw.XJitterWidth);
    
    
end

xlim([0 numel(bands_info) + 1])
xticks(1:numel(bands_info))
xticklabels({bands_info.name})
ylim([-1 15])
title('Average PLV difference (in %) for each channel')
ylabel('% Variation','Interpreter','none')
set(gca,'TickLabelInterpreter','none')



end


function plot_corr_channels(plv_eeg_expert_dataset,plv_ETL_dataset,bands_info, channels)

figure('WindowState', 'maximized');
hold on

% For each band
for iband = 1 : numel(bands_info)
    
    % Select the band data
    current_plv_global_eeg_expert = squeeze(plv_eeg_expert_dataset(:,:,iband,:));
    current_plv_global_ETL = squeeze(plv_ETL_dataset(:,:,iband,:));
    plv_vector_eeg_expert = reshape(current_plv_global_eeg_expert,[],size(current_plv_global_eeg_expert,3));
    plv_vector_ETL = reshape(current_plv_global_ETL,[],size(current_plv_global_ETL,3));
    
    % Get the upper matrix
    dummy = ones(size(current_plv_global_ETL(:,:,1)));
    upper_index = abs(triu(dummy,1)) > 0 ;
    upper_index = upper_index(:);
    plv_vector_eeg_expert =  plv_vector_eeg_expert(upper_index,:);
    plv_vector_ETL = plv_vector_ETL(upper_index,:);
    
    % Estimate correlation for each channel
    plv_vector_eeg_expert = plv_vector_eeg_expert';
    plv_vector_ETL = plv_vector_ETL';
    rho = nan(1,size(plv_vector_ETL,2));
    x_vector = iband * ones(1,size(plv_vector_ETL,2));
    for ichannel = 1 : size(plv_vector_ETL,2)
        
        % Remove nans
        current_channel_eeg_expert = plv_vector_eeg_expert(:,ichannel);
        current_channel_ETL = plv_vector_ETL(:,ichannel);
        nan_index = isnan(current_channel_eeg_expert) | isnan(current_channel_ETL);
        current_channel_eeg_expert = current_channel_eeg_expert(~nan_index);
        current_channel_ETL = current_channel_ETL(~nan_index);
        
        
        rho(ichannel) = corr(current_channel_eeg_expert,current_channel_ETL);
        
    end
    
    % Plot
    sw = swarmchart(x_vector,rho,'filled','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5);
    bx = boxchart(x_vector,rho,'BoxFaceColor',sw.CData ./ 1.2,'WhiskerLineColor',sw.CData ./ 1.2,...
        'MarkerStyle','none','BoxWidth',sw.XJitterWidth);
    
    
end

xlim([0 numel(bands_info) + 1])
xticks(1:numel(bands_info))
xticklabels({bands_info.name})
title('Average correlation of PLV for each channel')
ylabel('rho','Interpreter','none')
set(gca,'TickLabelInterpreter','none')



end


function plot_corr_plv_subjects(plv_eeg_expert_dataset,plv_ETL_dataset,bands_info, channels_eeg_expert,channels_ETL)

figure('WindowState', 'maximized');
hold on

% For each band
for iband = 1 : numel(bands_info)
    
    % Create the rho empty and the plot axis
    rho = nan(1,size(plv_eeg_expert_dataset,4));
    x_vector = iband * ones(1,size(rho,2));
    for irho = 1 :numel(rho)
        
        % Get the current plv matrix
        current_plv_eeg_expert = plv_eeg_expert_dataset(:,:,iband,irho);
        current_plv_ETL = plv_ETL_dataset(:,:,iband,irho);
        
        % Create an index of valid channels
        valid = channels_eeg_expert(irho).channels_included_index & channels_ETL(irho).channels_included_index;
        
        % Remove badchannels
        current_plv_eeg_expert = current_plv_eeg_expert(valid,valid);
        current_plv_ETL = current_plv_ETL(valid,valid);
        
        % Get the upper matrix in vector format
        dummy = ones(size(current_plv_ETL));
        upper_index = abs(triu(dummy,1)) > 0 ;
        current_plv_eeg_expert = current_plv_eeg_expert(upper_index);
        current_plv_ETL = current_plv_ETL(upper_index);
        
        % Correlation
        rho(irho) = corr(current_plv_eeg_expert,current_plv_ETL);
        
    end
    
    % Plot
    sw = swarmchart(x_vector,rho,'filled','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5);
    bx = boxchart(x_vector,rho,'BoxFaceColor',sw.CData ./ 1.2,'WhiskerLineColor',sw.CData ./ 1.2,...
        'MarkerStyle','none','BoxWidth',sw.XJitterWidth);
    
    
end

xlim([0 numel(bands_info) + 1])
xticks(1:numel(bands_info))
xticklabels({bands_info.name})
title('Average correlation of PLV matrix for each subject')
ylabel('rho','Interpreter','none')
set(gca,'TickLabelInterpreter','none')



end


