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

% Load the power
[plv_SRM_dataset, ~, channels_SRM] = read_plv_dataset(dataset,'SRM_database');
[plv_ETL_dataset, bands_info, channels_ETL] = read_plv_dataset(dataset,'ETL_database');

% Plot difference in channels
plot_diff_channels(plv_SRM_dataset,plv_ETL_dataset,bands_info)

% Plot correlation of each pair of channels
% plot_corr_channels(plv_SRM_dataset,plv_ETL_dataset,bands_info, channels_ETL)

% Plot correlation of PLV matrix for each subject
plot_corr_plv_subjects(plv_SRM_dataset,plv_ETL_dataset,bands_info, channels_SRM,channels_ETL)

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


function plot_diff_channels(plv_SRM_dataset,plv_ETL_dataset,bands_info)

figure('WindowState', 'maximized');
hold on

% For each band
for iband = 1 : numel(bands_info)
    
    % Select the band data
    current_plv_global_SRM = squeeze(plv_SRM_dataset(:,:,iband,:));
    current_plv_global_ETL = squeeze(plv_ETL_dataset(:,:,iband,:));
    
    % Get the difference of each channel and average all the differences
    plv_diff = current_plv_global_SRM - current_plv_global_ETL;
    plv_diff_avg = nanmean(plv_diff,3);
    
    % Symmetrical. Get only upper matrix
    plv_diff_vector = triu(ones(size(plv_diff_avg)),1) > 0 ;
    plv_diff_vector = plv_diff_avg(plv_diff_vector);
    x_vector = iband * ones(numel(plv_diff_vector),1);
    
    % Plot
    sw = swarmchart(x_vector,plv_diff_vector,'filled','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5);
    bx = boxchart(x_vector,plv_diff_vector,'BoxFaceColor',sw.CData ./ 1.2,'WhiskerLineColor',sw.CData ./ 1.2,...
        'MarkerStyle','none','BoxWidth',sw.XJitterWidth);
    
    
end

xlim([0 numel(bands_info) + 1])
xticks(1:numel(bands_info))
xticklabels({bands_info.name})
title('Average PLV difference for each channel')
ylim([-0.1 0.1])
ylabel('SRM_PLV - ETL_PLV','Interpreter','none')
set(gca,'TickLabelInterpreter','none')



end


function plot_corr_channels(plv_SRM_dataset,plv_ETL_dataset,bands_info, channels)

figure('WindowState', 'maximized');
hold on

% For each band
for iband = 1 : numel(bands_info)
    
    % Select the band data
    current_plv_global_SRM = squeeze(plv_SRM_dataset(:,:,iband,:));
    current_plv_global_ETL = squeeze(plv_ETL_dataset(:,:,iband,:));
    plv_vector_SRM = reshape(current_plv_global_SRM,[],size(current_plv_global_SRM,3));
    plv_vector_ETL = reshape(current_plv_global_ETL,[],size(current_plv_global_ETL,3));
    
    % Get the upper matrix
    dummy = ones(size(current_plv_global_ETL(:,:,1)));
    upper_index = abs(triu(dummy,1)) > 0 ;
    upper_index = upper_index(:);
    plv_vector_SRM =  plv_vector_SRM(upper_index,:);
    plv_vector_ETL = plv_vector_ETL(upper_index,:);
    
    % Estimate correlation for each channel
    plv_vector_SRM = plv_vector_SRM';
    plv_vector_ETL = plv_vector_ETL';
    rho = nan(1,size(plv_vector_ETL,2));
    x_vector = iband * ones(1,size(plv_vector_ETL,2));
    for ichannel = 1 : size(plv_vector_ETL,2)
        
        % Remove nans
        current_channel_SRM = plv_vector_SRM(:,ichannel);
        current_channel_ETL = plv_vector_ETL(:,ichannel);
        nan_index = isnan(current_channel_SRM) | isnan(current_channel_ETL);
        current_channel_SRM = current_channel_SRM(~nan_index);
        current_channel_ETL = current_channel_ETL(~nan_index);
        
        
        rho(ichannel) = corr(current_channel_SRM,current_channel_ETL);
        
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
ylim([0 1])
set(gca,'TickLabelInterpreter','none')



end


function plot_corr_plv_subjects(plv_SRM_dataset,plv_ETL_dataset,bands_info, channels_SRM,channels_ETL)

figure('WindowState', 'maximized');
hold on

% For each band
for iband = 1 : numel(bands_info)
    
    % Create the rho empty and the plot axis
    rho = nan(1,size(plv_SRM_dataset,4));
    x_vector = iband * ones(1,size(rho,2));
    for irho = 1 :numel(rho)
        
        % Get the current plv matrix
        current_plv_SRM = plv_SRM_dataset(:,:,iband,irho);
        current_plv_ETL = plv_ETL_dataset(:,:,iband,irho);
        
        % Create an index of valid channels
        valid = channels_SRM(irho).channels_included_index & channels_ETL(irho).channels_included_index;
        
        % Remove badchannels
        current_plv_SRM = current_plv_SRM(valid,valid);
        current_plv_ETL = current_plv_ETL(valid,valid);
        
        % Get the upper matrix in vector format
        dummy = ones(size(current_plv_ETL));
        upper_index = abs(triu(dummy,1)) > 0 ;
        current_plv_SRM = current_plv_SRM(upper_index);
        current_plv_ETL = current_plv_ETL(upper_index);
        
        % Correlation
        rho(irho) = corr(current_plv_SRM,current_plv_ETL);
        
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


