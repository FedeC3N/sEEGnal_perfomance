clear
clc
close all
restoredefaultpath

% Paths
config.path.dataset = '../../../../data/SRM_database/dataset';
config.path.stats = '../../../../data/SRM_database/stats';

% Create output path
if ~exist(config.path.stats), mkdir(config.path.stats), end

% Load the whole dataset and the stats
load(sprintf('%s/SRM_dataset.mat',config.path.dataset));
load(sprintf('%s/pow_stats.mat',config.path.stats));

% Load the power
[pow_SRM_dataset,~,channels_SRM] = read_pow_dataset(dataset,'SRM_database');
[pow_ETL_dataset,f,channels_ETL] = read_pow_dataset(dataset,'ETL_database');

%%%%%%%%%%%%%%%
% PLOTS
%%%%%%%%%%%%%%%

% Normalized pow spectrum for a specific band-area
% Prepare the pow_spectrum
to_plot = [];
to_plot.band_of_interest = {'broadband'};
to_plot.area_of_interest = {'whole_head'};
to_plot.areas_info = areas_info;
to_plot.bands_info = bands_info;

plot_pow_spectrum_norm_band_area(to_plot,pow_SRM_dataset,pow_ETL_dataset);

% Error bars in each area
% Prepare the pow_spectrum
to_plot = [];
to_plot.band_of_interest = {'broadband'};
to_plot.areas_info = areas_info;
to_plot.bands_info = bands_info;

plot_areas_errorbar(to_plot,pow_SRM_dataset,pow_ETL_dataset);

% Normalized pow spectrum for all areas in broadband
% Prepare the pow_spectrum
to_plot = [];
to_plot.areas_info = areas_info;
to_plot.bands_info = bands_info;

plot_pow_spectrum_norm_all_areas(to_plot,pow_SRM_dataset,pow_ETL_dataset);



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
    channels(icurrent).channels_included_index = ismember(pow.complete_channel_labels,pow.channels_included);
   
end

end


function plot_pow_spectrum_norm_band_area(to_plot,pow_SRM_dataset,pow_ETL_dataset)

% All channels
complete_channel_labels = {'Fp1';'AF7';'AF3';'F1';'F3';'F5';'F7';'FT7';'FC5';'FC3';'FC1';'C1';'C3';'C5';'T7';'TP7';'CP5';'CP3';'CP1';'P1';'P3';'P5';'P7';'P9';'PO7';'PO3';'O1';'Iz';'Oz';'POz';'Pz';'CPz';'Fpz';'Fp2';'AF8';'AF4';'AFz';'Fz';'F2';'F4';'F6';'F8';'FT8';'FC6';'FC4';'FC2';'FCz';'Cz';'C2';'C4';'C6';'T8';'TP8';'CP6';'CP4';'CP2';'P2';'P4';'P6';'P8';'P10';'PO8';'PO4';'O2'};

% Find electrodes of interest
areas_index = find(ismember({to_plot.areas_info.name},to_plot.area_of_interest));
areas_of_interest = to_plot.areas_info(areas_index).channel;
channels_index = ismember(complete_channel_labels,areas_of_interest);

% Find band of interest
band_index = find(ismember({to_plot.bands_info.name},to_plot.band_of_interest));
f_index = to_plot.bands_info(band_index).f_limits_index;

% f of interest
f_of_interest = to_plot.bands_info(band_index).f_original;
f_of_interest = f_of_interest(f_index);

% Update pow matrix to plot
pow_SRM_dataset = pow_SRM_dataset(channels_index,f_index,:);
pow_ETL_dataset = pow_ETL_dataset(channels_index,f_index,:);

% Average across electrodes
SRM_average = squeeze(nanmean(nanmean(pow_SRM_dataset,1),3));
ETL_average = squeeze(nanmean(nanmean(pow_ETL_dataset,1),3));

% Estimate mean error
SRM_std_error = nanstd(squeeze(nanmean(pow_SRM_dataset,1)),1,2)/sqrt(size(pow_SRM_dataset,3));
ETL_std_error = nanstd(squeeze(nanmean(pow_ETL_dataset,1)),1,2)/sqrt(size(pow_ETL_dataset,3));

% Color
plot_color_SRM = [255,179,179]/255;
plot_color_ETL = [179,179,255]/255;

% Plot
figure
patch([f_of_interest,fliplr(f_of_interest)], ...
    [(SRM_average - SRM_std_error'),fliplr((SRM_average + SRM_std_error'))],...
    plot_color_SRM,'FaceAlpha',0.2,'HandleVisibility','off')
hold on
plot(f_of_interest,SRM_average,'r','LineWidth',2)
patch([f_of_interest,fliplr(f_of_interest)], ...
    [(ETL_average - ETL_std_error'),fliplr((ETL_average + ETL_std_error'))],...
    plot_color_ETL,'FaceAlpha',0.2,'HandleVisibility','off')
plot(f_of_interest,ETL_average,'b','LineWidth',2)

% Enhance the plot
legend({'SRM database', 'ETL database'});
title(sprintf('Power in %s and %s area', to_plot.band_of_interest{1}, to_plot.area_of_interest{1}),...
    'Interpreter','none')

end


function plot_areas_errorbar(to_plot,pow_SRM_dataset,pow_ETL_dataset)

% All channels
complete_channel_labels = {'Fp1';'AF7';'AF3';'F1';'F3';'F5';'F7';'FT7';'FC5';'FC3';'FC1';'C1';'C3';'C5';'T7';'TP7';'CP5';'CP3';'CP1';'P1';'P3';'P5';'P7';'P9';'PO7';'PO3';'O1';'Iz';'Oz';'POz';'Pz';'CPz';'Fpz';'Fp2';'AF8';'AF4';'AFz';'Fz';'F2';'F4';'F6';'F8';'FT8';'FC6';'FC4';'FC2';'FCz';'Cz';'C2';'C4';'C6';'T8';'TP8';'CP6';'CP4';'CP2';'P2';'P4';'P6';'P8';'P10';'PO8';'PO4';'O2'};

avg_pow = nan(2,numel(to_plot.areas_info));
std_error_pow = nan(2,numel(to_plot.areas_info));

for iarea = 1 : numel(to_plot.areas_info)
    
    % Find electrodes of interest
    areas_of_interest = to_plot.areas_info(iarea).channel;
    channels_index = ismember(complete_channel_labels,areas_of_interest);
    
    % Get pow matrix to plot
    current_pow_SRM_dataset = pow_SRM_dataset(channels_index,:,:);
    current_pow_ETL_dataset = pow_ETL_dataset(channels_index,:,:);
    
    % Average across electrodes and subjects
    current_pow_SRM_dataset = squeeze(nanmean(current_pow_SRM_dataset,1));
    
    current_pow_SRM_dataset = nanmean(squeeze(nanmean(current_pow_SRM_dataset,1)),1);
    current_pow_ETL_dataset = squeeze(nanmean(nanmean(current_pow_ETL_dataset,1),2));
    current_pow_SRM_dataset = current_pow_SRM_dataset';
    current_pow_ETL_dataset = current_pow_ETL_dataset';
    
    % Estimate the average and mean error
    SRM_std_error = nanstd(current_pow_SRM_dataset)/sqrt(numel(current_pow_SRM_dataset));
    ETL_std_error = nanstd(current_pow_ETL_dataset)/sqrt(numel(current_pow_ETL_dataset));

    % Save
    avg_pow(1,iarea) = nanmean(current_pow_SRM_dataset);
    avg_pow(2,iarea) = nanmean(current_pow_ETL_dataset);
    std_error_pow(1,iarea) = SRM_std_error;
    std_error_pow(2,iarea) = ETL_std_error;
    
    
end

% Reshape avg_pow for plotting purposes
avg_pow = avg_pow';
std_error_pow = std_error_pow';

% First the bar plot
figure
hb = bar(avg_pow);
hold on


% Plot the error bars
for igroup = 1:size(avg_pow,2)
    % get x positions per group
    xpos = hb(igroup).XData + hb(igroup).XOffset;
    % draw errorbar
    errorbar(xpos, avg_pow(:,igroup), std_error_pow(:,igroup), 'LineStyle', 'none', ... 
        'Color', 'k', 'LineWidth', 1);
end

% Enhance the plot
legend({'SRM database', 'ETL database'});
xticklabels({to_plot.areas_info.name})
set(gca,'TickLabelInterpreter','none')
title('Average power per area across')



end


function plot_pow_spectrum_norm_all_areas(to_plot,pow_SRM_dataset,pow_ETL_dataset)

% All channels
complete_channel_labels = {'Fp1';'AF7';'AF3';'F1';'F3';'F5';'F7';'FT7';'FC5';'FC3';'FC1';'C1';'C3';'C5';'T7';'TP7';'CP5';'CP3';'CP1';'P1';'P3';'P5';'P7';'P9';'PO7';'PO3';'O1';'Iz';'Oz';'POz';'Pz';'CPz';'Fpz';'Fp2';'AF8';'AF4';'AFz';'Fz';'F2';'F4';'F6';'F8';'FT8';'FC6';'FC4';'FC2';'FCz';'Cz';'C2';'C4';'C6';'T8';'TP8';'CP6';'CP4';'CP2';'P2';'P4';'P6';'P8';'P10';'PO8';'PO4';'O2'};

% New figure to plot
figure

% Go through each area
for iarea = 1 : numel(to_plot.areas_info)
    
    % Find electrodes of interest
    areas_of_interest = to_plot.areas_info(iarea).channel;
    channels_index = ismember(complete_channel_labels,areas_of_interest);
    
    % Select broadban frequencies
    f_index = to_plot.bands_info(6).f_limits_index;
    f_of_interest = to_plot.bands_info(6).f_original;
    f_of_interest = f_of_interest(f_index);
    
    % Update pow matrix to plot
    current_pow_SRM_dataset = pow_SRM_dataset(channels_index,f_index,:);
    current_pow_ETL_dataset = pow_ETL_dataset(channels_index,f_index,:);
    
    % Average across electrodes
    SRM_average = squeeze(nanmean(nanmean(current_pow_SRM_dataset,1),3));
    ETL_average = squeeze(nanmean(nanmean(current_pow_ETL_dataset,1),3));
    
    % Estimate mean error
    SRM_std_error = nanstd(squeeze(nanmean(current_pow_SRM_dataset,1)),1,2)/sqrt(size(current_pow_SRM_dataset,3));
    ETL_std_error = nanstd(squeeze(nanmean(current_pow_ETL_dataset,1)),1,2)/sqrt(size(current_pow_ETL_dataset,3));
    
    % Color
    plot_color_SRM = [255,179,179]/255;
    plot_color_ETL = [179,179,255]/255;
    
    % Plot
    subplot(2,4,iarea)
    patch([f_of_interest,fliplr(f_of_interest)], ...
        [(SRM_average - SRM_std_error'),fliplr((SRM_average + SRM_std_error'))],...
        plot_color_SRM,'FaceAlpha',0.2,'HandleVisibility','off')
    hold on
    plot(f_of_interest,SRM_average,'r','LineWidth',2)
    patch([f_of_interest,fliplr(f_of_interest)], ...
        [(ETL_average - ETL_std_error'),fliplr((ETL_average + ETL_std_error'))],...
        plot_color_ETL,'FaceAlpha',0.2,'HandleVisibility','off')
    plot(f_of_interest,ETL_average,'b','LineWidth',2)
    ylim([0 0.0008])
    
    % Enhance the plot
    legend({'SRM database', 'ETL database'});
    title(sprintf('Power in BroadBand and %s area',  to_plot.areas_info(iarea).name),...
        'Interpreter','none')
    
end

end
