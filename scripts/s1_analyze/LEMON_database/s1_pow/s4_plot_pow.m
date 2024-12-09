clear
clc
close all
restoredefaultpath

% Paths
config.path.clean_data = '../../../../databases/LEMON_database/derivatives';
config.path.results = '../../../../results/pow';

% Load the results
load(sprintf('%s/pow_results.mat',config.path.results));

% Get the different testers
testers = {'lemon','sEEGnal'};

% Channels
config.complete_channel_labels = {'Fp1', 'Fp2', 'F7', 'F3', 'Fz', 'F4', 'F8', 'FC5', 'FC1', 'FC2', 'FC6', 'T7', 'C3', 'Cz', 'C4', 'T8', 'CP5', 'CP1', 'CP2', 'CP6', 'AFz', 'P7', 'P3', 'Pz', 'P4', 'P8', 'PO9', 'O1', 'Oz', 'O2', 'PO10', 'AF7', 'AF3', 'AF4', 'AF8', 'F5', 'F1', 'F2', 'F6', 'FT7', 'FC3', 'FC4', 'FT8', 'C5', 'C1', 'C2', 'C6', 'TP7', 'CP3', 'CPz', 'CP4', 'TP8', 'P5', 'P1', 'P2', 'P6', 'PO7', 'PO3', 'POz', 'PO4', 'PO8'};

% Read the power spectrum
[pow_lemon_dataset,~,channels_lemon_included] = read_pow_dataset(config,'lemon');
[pow_sEEGnal_dataset,f,channels_sEEGnal_included] = read_pow_dataset(config,'sEEGnal');

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

plot_pow_spectrum_norm_band_area(config,to_plot,pow_lemon_dataset,pow_sEEGnal_dataset);

% Error bars in each area
% Prepare the pow_spectrum
to_plot = [];
to_plot.band_of_interest = {'broadband'};
to_plot.areas_info = areas_info;
to_plot.bands_info = bands_info;

plot_areas_errorbar(config,to_plot,pow_lemon_dataset,pow_sEEGnal_dataset);

% Normalized pow spectrum for all areas in broadband
% Prepare the pow_spectrum
to_plot = [];
to_plot.areas_info = areas_info;
to_plot.bands_info = bands_info;

plot_pow_spectrum_norm_all_areas(config,to_plot,pow_lemon_dataset,pow_sEEGnal_dataset);



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

function plot_pow_spectrum_norm_band_area(config,to_plot,pow_lemon_dataset,pow_sEEGnal_dataset)

% All channels
complete_channel_labels = config.complete_channel_labels;

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
pow_lemon_dataset = pow_lemon_dataset(channels_index,f_index,:);
pow_sEEGnal_dataset = pow_sEEGnal_dataset(channels_index,f_index,:);

% Average across electrodes
lemon_average = squeeze(nanmean(nanmean(pow_lemon_dataset,1),3));
sEEGnal_average = squeeze(nanmean(nanmean(pow_sEEGnal_dataset,1),3));

% Estimate mean error
lemon_std_error = nanstd(squeeze(nanmean(pow_lemon_dataset,1)),1,2)/sqrt(size(pow_lemon_dataset,3));
sEEGnal_std_error = nanstd(squeeze(nanmean(pow_sEEGnal_dataset,1)),1,2)/sqrt(size(pow_sEEGnal_dataset,3));

% Color
plot_color_lemon = [255,179,179]/255;
plot_color_sEEGnal = [179,179,255]/255;

% Plot
figure('WindowState','maximized')
patch([f_of_interest,fliplr(f_of_interest)], ...
    [(lemon_average - lemon_std_error'),fliplr((lemon_average + lemon_std_error'))],...
    plot_color_lemon,'FaceAlpha',0.2,'HandleVisibility','off')
hold on
plot(f_of_interest,lemon_average,'r','LineWidth',2)
patch([f_of_interest,fliplr(f_of_interest)], ...
    [(sEEGnal_average - sEEGnal_std_error'),fliplr((sEEGnal_average + sEEGnal_std_error'))],...
    plot_color_sEEGnal,'FaceAlpha',0.2,'HandleVisibility','off')
plot(f_of_interest,sEEGnal_average,'b','LineWidth',2)

% Enhance the plot
legend({'Lemon database', 'sEEGnal database'});
title(sprintf('Power in %s and %s area', to_plot.band_of_interest{1}, to_plot.area_of_interest{1}),...
    'Interpreter','none')

end

function plot_areas_errorbar(config,to_plot,pow_lemon_dataset,pow_sEEGnal_dataset)

% All channels
complete_channel_labels = config.complete_channel_labels;

avg_pow = nan(2,numel(to_plot.areas_info));
std_error_pow = nan(2,numel(to_plot.areas_info));

for iarea = 1 : numel(to_plot.areas_info)
    
    % Find electrodes of interest
    areas_of_interest = to_plot.areas_info(iarea).channel;
    channels_index = ismember(complete_channel_labels,areas_of_interest);
    
    % Get pow matrix to plot
    current_pow_lemon_dataset = pow_lemon_dataset(channels_index,:,:);
    current_pow_sEEGnal_dataset = pow_sEEGnal_dataset(channels_index,:,:);
    
    % Average across electrodes and subjects
    current_pow_lemon_dataset = squeeze(nanmean(current_pow_lemon_dataset,1));
    
    current_pow_lemon_dataset = nanmean(squeeze(nanmean(current_pow_lemon_dataset,1)),1);
    current_pow_sEEGnal_dataset = squeeze(nanmean(nanmean(current_pow_sEEGnal_dataset,1),2));
    current_pow_lemon_dataset = current_pow_lemon_dataset';
    current_pow_sEEGnal_dataset = current_pow_sEEGnal_dataset';
    
    % Estimate the average and mean error
    lemon_std_error = nanstd(current_pow_lemon_dataset)/sqrt(numel(current_pow_lemon_dataset));
    sEEGnal_std_error = nanstd(current_pow_sEEGnal_dataset)/sqrt(numel(current_pow_sEEGnal_dataset));

    % Save
    avg_pow(1,iarea) = nanmean(current_pow_lemon_dataset);
    avg_pow(2,iarea) = nanmean(current_pow_sEEGnal_dataset);
    std_error_pow(1,iarea) = lemon_std_error;
    std_error_pow(2,iarea) = sEEGnal_std_error;
    
    
end

% Reshape avg_pow for plotting purposes
avg_pow = avg_pow';
std_error_pow = std_error_pow';

% First the bar plot
figure('WindowState','maximized')
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
legend({'Lemon database', 'sEEGnal database'});
xticklabels({to_plot.areas_info.name})
set(gca,'TickLabelInterpreter','none')
title('Average power per area across')



end

function plot_pow_spectrum_norm_all_areas(config,to_plot,pow_lemon_dataset,pow_sEEGnal_dataset)

% All channels
complete_channel_labels = config.complete_channel_labels;

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
    current_pow_lemon_dataset = pow_lemon_dataset(channels_index,f_index,:);
    current_pow_sEEGnal_dataset = pow_sEEGnal_dataset(channels_index,f_index,:);
    
    % Average across electrodes
    lemon_average = squeeze(nanmean(nanmean(current_pow_lemon_dataset,1),3));
    sEEGnal_average = squeeze(nanmean(nanmean(current_pow_sEEGnal_dataset,1),3));
    
    % Estimate mean error
    lemon_std_error = nanstd(squeeze(nanmean(current_pow_lemon_dataset,1)),1,2)/sqrt(size(current_pow_lemon_dataset,3));
    sEEGnal_std_error = nanstd(squeeze(nanmean(current_pow_sEEGnal_dataset,1)),1,2)/sqrt(size(current_pow_sEEGnal_dataset,3));
    
    % Color
    plot_color_lemon = [255,179,179]/255;
    plot_color_sEEGnal = [179,179,255]/255;
    
    % Plot
    figure('WindowState','maximized')
    patch([f_of_interest,fliplr(f_of_interest)], ...
        [(lemon_average - lemon_std_error'),fliplr((lemon_average + lemon_std_error'))],...
        plot_color_lemon,'FaceAlpha',0.2,'HandleVisibility','off')
    hold on
    plot(f_of_interest,lemon_average,'r','LineWidth',2)
    patch([f_of_interest,fliplr(f_of_interest)], ...
        [(sEEGnal_average - sEEGnal_std_error'),fliplr((sEEGnal_average + sEEGnal_std_error'))],...
        plot_color_sEEGnal,'FaceAlpha',0.2,'HandleVisibility','off')
    plot(f_of_interest,sEEGnal_average,'b','LineWidth',2)
    
    % Enhance the plot
    legend({'Lemon database', 'sEEGnal database'});
    title(sprintf('Power in BroadBand and %s area',  to_plot.areas_info(iarea).name),...
        'Interpreter','none')
    
end

end
