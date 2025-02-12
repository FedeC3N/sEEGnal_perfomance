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

% Define the frequency bands
config.bands = {'delta', 'theta','alpha','beta','gamma', 'broadband'};

% Define measures
% config.measures = {'NRMSE', 'rho', 'tstat'};
config.measures = {'NRMSE', 'rho'};

% Read the power spectrum
[pow_lemon_dataset,~,channels_lemon_included] = read_pow_dataset(config,'lemon');
[pow_sEEGnal_dataset,f,channels_sEEGnal_included] = read_pow_dataset(config,'sEEGnal');

%%%%%%%%%%%%%%%
% PLOTS
%%%%%%%%%%%%%%%

% Topoplots of the stats
% to_plot = [];
% to_plot.bands_info = bands_info;
% 
% plot_stats_in_head(config,to_plot,stats,pow_lemon_dataset,pow_sEEGnal_dataset)

% Normalized pow spectrum for all areas in broadband
% Prepare the pow_spectrum
to_plot = [];
to_plot.areas_info = areas_info;
to_plot.bands_info = bands_info;

plot_pow_spectrum_norm_all_areas(config,to_plot,stats,pow_lemon_dataset,pow_sEEGnal_dataset);

% Plot lines between lemon and sEEGnal
% plot_lines_linking_values(pow_lemon_dataset,pow_sEEGnal_dataset,bands_info)

% Normalized pow spectrum for a specific band-area
% Prepare the pow_spectrum
% to_plot = [];
% to_plot.band_of_interest = {'broadband'};
% to_plot.area_of_interest = {'whole_head'};
% to_plot.areas_info = areas_info;
% to_plot.bands_info = bands_info;
% 
% plot_pow_spectrum_norm_band_area(config,to_plot,pow_lemon_dataset,pow_sEEGnal_dataset);

% Error bars in each area
% Prepare the pow_spectrum
% to_plot = [];
% to_plot.band_of_interest = {'broadband'};
% to_plot.areas_info = areas_info;
% to_plot.bands_info = bands_info;
% 
% plot_areas_errorbar(config,to_plot,pow_lemon_dataset,pow_sEEGnal_dataset);





% Functions
function [pow_dataset_norm,f,channels] = read_pow_dataset(config,dataset_name)

% Load the datset
dataset_path = sprintf('%s/%s/%s_dataset.mat',config.path.clean_data,...
    dataset_name,dataset_name);
dummy = load(dataset_path);

for icurrent = 1 : numel(dummy.dataset)
    
    % Load pow
    current_dataset = dummy.dataset(icurrent);
    pow_file = sprintf('../../../../%s/%s',current_dataset.pow.path,...
        current_dataset.pow.file);
    pow = load(pow_file);
    
    % Save f
    f = pow.f;
    
    % Normalize
    current_pow = nanmean(pow.pow_spectrum,3);
    scaling_factor = 1./(nansum(current_pow,2));
    scaling_factor = repmat(scaling_factor,[1 size(current_pow,2)]);
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


function plot_stats_in_head(config,to_plot,stats,pow_lemon_dataset,pow_sEEGnal_dataset)

% Draw for each band and measure
for iband = 1 : numel(config.bands)
    
    current_band = config.bands{iband};
    
    % Define figure
    fig = figure('WindowState', 'maximized');
    
    for imeasure = 1 : numel(config.measures)
        
        current_measure = config.measures{imeasure};
        
        % Scatter the head
        ax1 = subplot(2,2,imeasure);
        [pos_elec, size_elec] = draw_head(config);
        
        % Scatter the values of interest
        color_elec = nanmean(stats.(current_band).(current_measure),2);
        scatter(pos_elec(:,1),pos_elec(:,2),size_elec,color_elec,'filled')
        c = colorbar;
        
        % Set limits to the colorbar
%         if min(color_elec) < 0
%             caxis([min(color_elec), max(color_elec)]);
%         else
%             caxis([0, max(color_elec)]);
%         end
        
        % Set the colormap
        gradient_1_neg = linspace(1,1,2000);
        gradient_2_neg = linspace(1,0.4,2000);
        gradient_3_neg = linspace(1,0.4,2000);
        gradient_1_pos = linspace(0.4,1,2000);
        gradient_2_pos = linspace(0.4,1,2000);
        gradient_3_pos = linspace(1,1,2000);
        clrmap = [[gradient_1_pos gradient_1_neg]',...
            [gradient_2_pos gradient_2_neg]',...
            [gradient_3_pos gradient_3_neg]'];
        if min(color_elec) > 0
            clrmap = clrmap(2000:end,:);
        end
        colormap(ax1,clrmap)
        
        
        
        
        % Details
        title(sprintf('Band %s - Measure %s', current_band, current_measure))
        
        
    end
 
    
    % Select broadban frequencies
    f_index = to_plot.bands_info(iband).f_limits_index;
    f_of_interest = to_plot.bands_info(iband).f_original;
    f_of_interest = f_of_interest(f_index);
    
    % Update pow matrix to plot
    current_pow_lemon_dataset = pow_lemon_dataset(:,f_index,:);
    current_pow_sEEGnal_dataset = pow_sEEGnal_dataset(:,f_index,:);
    
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
    ax1 = subplot(2,2,4);
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
    title(sprintf('Power in band %s and the whole head',current_band),'Interpreter','none')
    
    
end

end


function plot_lines_linking_values(pow_lemon_dataset,pow_sEEGnal_dataset,bands_info)

figure('WindowState', 'maximized');
hold on

% Color
plot_colors = [179,179,255;...
    255,144,144;...
    20,144,144;...
    20,70,144;...
    20,70,70;...
    255,20,144]/255;

for iband = 1 : numel(bands_info)
    
    
    
    % Get the current band
    current_f = bands_info(iband).f_limits_index;
    
    % Estimate the plots
    lemon_average = nanmean(nanmean(pow_lemon_dataset(:,current_f,:),3),2);
    lemon_average = normalize(lemon_average,'range');
    lemon_average = nanmean(lemon_average);
    sEEGnal_average = nanmean(nanmean(pow_sEEGnal_dataset(:,current_f,:),3),2);
    sEEGnal_average = normalize(sEEGnal_average,'range');
    sEEGnal_average = nanmean(sEEGnal_average);
    lemon_std_error = nanmean(nanmean(pow_lemon_dataset,3),2);
    lemon_std_error = normalize(lemon_std_error,'range');
    lemon_std_error = lemon_std_error/sqrt(size(pow_lemon_dataset,1));
    lemon_std_error = nanmean(lemon_std_error);
    sEEGnal_std_error = nanmean(nanmean(pow_sEEGnal_dataset,3),2);
    sEEGnal_std_error = normalize(sEEGnal_std_error,'range');
    sEEGnal_std_error = sEEGnal_std_error/sqrt(size(pow_sEEGnal_dataset,1));
    sEEGnal_std_error = nanmean(sEEGnal_std_error);    
    
%     patch([[1,2],[2,1]], ...
%         [[lemon_average + lemon_std_error, sEEGnal_average + sEEGnal_std_error],...
%         [sEEGnal_average - sEEGnal_std_error, lemon_average - lemon_std_error]],...
%         plot_colors(iband,:),'FaceAlpha',0.2,'HandleVisibility','off')
    plot([1 2],[lemon_average, sEEGnal_average],'Color',plot_colors(iband,:),'LineWidth',2)    
          
end

% Plot characteristics
xlim([0.5 2.5])
ylim([0 1])
legend(bands_info.name)


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


function plot_pow_spectrum_norm_all_areas(config,to_plot,stats,pow_lemon_dataset,pow_sEEGnal_dataset)

% All channels
complete_channel_labels = config.complete_channel_labels;

% Go through each area
for iarea = 1 : numel(to_plot.areas_info)
    
    % Find electrodes of interest
    areas_of_interest = to_plot.areas_info(iarea).channel;
    channels_index = ismember(complete_channel_labels,areas_of_interest);
    
    % Define figure
    fig = figure('WindowState', 'maximized');
    
    for imeasure = 1 : numel(config.measures)
        
        current_measure = config.measures{imeasure};
        
        % Scatter the head
        ax1 = subplot(1,2,imeasure);
        [pos_elec, size_elec] = draw_head(config);
        color_elec = nanmean(stats.broadband.(current_measure),2);
        
        % Keep the channels of interest
        pos_elec = pos_elec(channels_index,:);
        size_elec = size_elec(channels_index);
        color_elec = color_elec(channels_index,:);
        
        % Scatter the values of interest
        scatter(pos_elec(:,1),pos_elec(:,2),size_elec,color_elec,'filled')
        c = colorbar;
        
        % Set the colormap
        gradient_1_neg = linspace(1,1,2000);
        gradient_2_neg = linspace(1,0.4,2000);
        gradient_3_neg = linspace(1,0.4,2000);
        gradient_1_pos = linspace(0.4,1,2000);
        gradient_2_pos = linspace(0.4,1,2000);
        gradient_3_pos = linspace(1,1,2000);
        clrmap = [[gradient_1_pos gradient_1_neg]',...
            [gradient_2_pos gradient_2_neg]',...
            [gradient_3_pos gradient_3_neg]'];
        if min(color_elec) > 0
            clrmap = clrmap(2000:end,:);
        end
        colormap(ax1,clrmap)
    
        % Details
        title(sprintf('Broadband - Measure %s', current_measure))
             
    end
        
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

    hold on
    plot(f_of_interest,nanmean(current_pow_lemon_dataset,3),...
        'Color',[1 0.6 0.6],'LineWidth',0.5)
    plot(f_of_interest,nanmean(current_pow_sEEGnal_dataset,3),...
        'Color',[0.6 0.6 1],'LineWidth',0.5)
%     patch([f_of_interest,fliplr(f_of_interest)], ...
%         [(lemon_average - lemon_std_error'),fliplr((lemon_average + lemon_std_error'))],...
%         plot_color_lemon,'FaceAlpha',0.2,'HandleVisibility','off')
%     patch([f_of_interest,fliplr(f_of_interest)], ...
%         [(sEEGnal_average - sEEGnal_std_error'),fliplr((sEEGnal_average + sEEGnal_std_error'))],...
%         plot_color_sEEGnal,'FaceAlpha',0.2,'HandleVisibility','off')
    plot(f_of_interest,lemon_average,'r','LineWidth',2)
    plot(f_of_interest,sEEGnal_average,'b','LineWidth',2)
    
    % Enhance the plot
    legend({'Lemon database', 'sEEGnal database'});
    title(sprintf('Power in BroadBand and %s area',  to_plot.areas_info(iarea).name),...
        'Interpreter','none')
    
end

end


function [pos_elec,size_elec] = draw_head(config)

hold on;
axis equal off

% Head (big circle)
theta = linspace(0, 2*pi, 100);
r_head = 0.5; % radius of the head
x_head = r_head * cos(theta);
y_head = r_head * sin(theta);

% Ears (two small circles)
r_ear = 0.1; % radius of the ears
x_ear = r_ear * cos(theta);
y_ear = r_ear * sin(theta);

% Left ear
fill(x_ear - r_head/1.2, y_ear, [0.9 0.9 0.9]); % Slightly darker gray

% Right ear
fill(x_ear + r_head/1.2, y_ear, [0.9 0.9 0.9]);

% Nose (triangle)
x_nose = [-0.07, 0.07, 0]; % X-coordinates of the triangle
y_nose = [0.49, 0.49, 0.55]; % Y-coordinates of the triangle
fill(x_nose, y_nose, [0.9 0.9 0.9]); % Darker gray for nose

% Plot the head at the end
fill(x_head, y_head, [0.9 0.9 0.9]); % Light gray color

% Plot the sensors
% Save memory for the positions and colors
pos_elec = nan(numel(config.complete_channel_labels),2);
size_elec = 50*ones(numel(config.complete_channel_labels),1);

% Read the channels position file
lines = readlines('../shared/head_layouts/elec1005.lay');

% Go through each line
for iline = 1 : size(lines,1)
    
    % Get the channel position in the original "complete_channel_labels"
    % for consistency
    current_line = strsplit(lines(iline),' ');
    
    if size(current_line,2) > 1
        current_channel = current_line(6);
        current_channel_index = ismember(config.complete_channel_labels,current_channel);
        
        % If present, save the position and color
        if sum(current_channel_index) == 1
            
            % Position
            pos_elec(current_channel_index,1) = current_line(2);
            pos_elec(current_channel_index,2) = current_line(3);
            
        end
        
    end
    
end

scatter(pos_elec(:,1),pos_elec(:,2),size_elec,'MarkerEdgeColor',[0 0 0])
pos_labels = [pos_elec(:,1)-0.015, pos_elec(:,2)+0.03];
text(pos_labels(:,1),pos_labels(:,2),config.complete_channel_labels)

end