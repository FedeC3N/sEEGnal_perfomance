clear
clc
close all
restoredefaultpath

% Paths
config.path.clean_data = '../../../../databases/AI_Mind_database/derivatives';
config.path.results = '../../../../results/AI_Mind_database/pow';
config.path.figures = '../../../../docs/manuscript/figures/AI_Mind_database/pow_results';

if ~exist(config.path.figures), mkdir(config.path.figures), end

% Load the results
load(sprintf('%s/pow_results.mat',config.path.results));

% Get the different testers (except sEEGnal)
testers = dir(sprintf('%s/*',config.path.clean_data));
testers = testers(3:end-1);
testers = {testers.name};
config.testers = testers;

% Channels
config.complete_channel_labels = {'Fp1', 'Fpz', 'Fp2', 'F7', 'F3', 'Fz', 'F4', 'F8', 'FC5', 'FC1', 'FC2', 'FC6',...
    'M1', 'T7', 'C3', 'Cz', 'C4', 'T8', 'M2', 'CP5', 'CP1', 'CP2', 'CP6', 'P7', 'P3',...
    'Pz', 'P4', 'P8', 'POz', 'O1', 'O2', 'AF7', 'AF3', 'AF4', 'AF8', 'F5', 'F1', 'F2',...
    'F6', 'FC3', 'FCz', 'FC4', 'C5', 'C1', 'C2', 'C6', 'CP3', 'CP4', 'P5', 'P1', 'P2',...
    'P6', 'F9', 'PO3', 'PO4', 'F10', 'FT7', 'FT8', 'TP7', 'TP8', 'PO7', 'PO8', 'FT9',...
    'FT10', 'TPP9h', 'TPP10h', 'PO9', 'PO10', 'P9', 'P10', 'AFF1', 'AFz', 'AFF2', 'FFC5h',...
    'FFC3h', 'FFC4h', 'FFC6h', 'FCC5h', 'FCC3h', 'FCC4h', 'FCC6h', 'CCP5h', 'CCP3h', 'CCP4h',...
    'CCP6h', 'CPP5h', 'CPP3h', 'CPP4h', 'CPP6h', 'PPO1', 'PPO2', 'I1', 'Iz', 'I2', 'AFp3h', 'AFp4h',...
    'AFF5h', 'AFF6h', 'FFT7h', 'FFC1h', 'FFC2h', 'FFT8h', 'FTT9h', 'FTT7h', 'FCC1h', 'FCC2h', 'FTT8h',...
    'FTT10h', 'TTP7h', 'CCP1h', 'CCP2h', 'TTP8h', 'TPP7h', 'CPP1h', 'CPP2h', 'TPP8h', 'PPO9h', 'PPO5h',...
    'PPO6h', 'PPO10h', 'POO9h', 'POO3h', 'POO4h', 'POO10h', 'OI1h', 'OI2h'};

% Areas
areas_info = struct('name',{'frontal','temporal_l','temporal_r','parietal_l',...
    'parietal_r','occipital','whole_head'},'channel',[]);
areas_info(1).channel = {'Fp1', 'Fp2', 'F7', 'F3', 'Fz', 'F4', 'F8', 'AF7', 'AF3', 'AF4', 'AF8', 'F5', 'F1', 'F2', 'F6'};
areas_info(2).channel = {'T7', 'FT7', 'TP7'};
areas_info(3).channel = {'T8', 'FT8', 'TP8'};
areas_info(4).channel = {'CP5', 'CP1', 'P7', 'P3', 'TP7', 'CP3', 'P5', 'P1'};
areas_info(5).channel = {'CP2', 'CP6', 'P4', 'P8', 'CP4', 'TP8', 'P2', 'P6'};
areas_info(6).channel = {'PO9', 'O1', 'Oz', 'O2', 'PO10','PO7', 'PO3', 'POz', 'PO4', 'PO8'};
areas_info(7).channel = complete_channel_labels;

% Define the frequency bands
config.bands = {'delta', 'theta','alpha','beta','gamma', 'broadband'};

% Define measures
% config.measures = {'NRMSE', 'rho', 'tstat'};
config.measures = {'NRMSE', 'rho'};

% Read the power spectrum
[pow_human_dataset] = read_pow_dataset(config,testers);
[pow_sEEGnal_dataset] = read_pow_dataset(config,{'sEEGnal'});

%%%%%%%%%%%%%%%
% PLOTS
%%%%%%%%%%%%%%%
% Normalized pow spectrum for all areas in broadband
% Prepare the pow_spectrum
to_plot = [];
to_plot.areas_info = areas_info;
to_plot.bands_info = bands_info;

plot_pow_spectrum_norm_all_areas(config,to_plot,stats,pow_human_dataset,pow_sEEGnal_dataset);

% Functions
function [pow_dataset_norm] = read_pow_dataset(config,dataset_name)

for itester = 1 : numel(dataset_name)
    
    % Load the datset
    dataset_path = sprintf('%s/%s/%s_dataset.mat',config.path.clean_data,...
        dataset_name{itester},dataset_name{itester});
    dummy = load(dataset_path);
    
    for icurrent = 1 : numel(dummy.dataset)
        
        % Load pow
        current_dataset = dummy.dataset(icurrent);
        pow_file = sprintf('../../../../%s/%s.mat',current_dataset.pow.path,...
            current_dataset.pow.file);
        
        if ~exist(pow_file)
            fake_pow = nan(size(current_pow));
            pow_dataset_norm(:,:,icurrent,itester) = fake_pow;
            continue
        end
        
        pow = load(pow_file);
        
        % Save f
        f = pow.f;
        
        % Normalize
        current_pow = nanmean(pow.pow_spectrum,3);
        scaling_factor = 1./(nansum(current_pow,2));
        scaling_factor = repmat(scaling_factor,[1 size(current_pow,2)]);
        current_pow_norm = scaling_factor .* current_pow;
        
        % Add to the all matrix
        if icurrent == 1 & itester == 1
            pow_dataset_norm = nan(numel(pow.complete_channel_labels),numel(f),numel(dummy.dataset),numel(dataset_name));
            channels = struct('channels_included',[],...
                'channels_included_index',[]);
        end
        pow_dataset_norm(1:size(current_pow_norm,1),:,icurrent,itester) = current_pow_norm;
        
    end
end

end


function plot_pow_spectrum_norm_all_areas(config,to_plot,stats,pow_human_dataset,pow_sEEGnal_dataset)

% All channels
complete_channel_labels = config.complete_channel_labels;

% Go through each area
for iarea = 7 : numel(to_plot.areas_info)
    
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
        current_stats = nanmean(stats.(current_measure)(:,:,6,:),4);
        color_elec = nanmean(current_stats,2);
        
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
    
    % Save the figure
    outfile = sprintf('%s/%s_head_measures.svg',config.path.figures,...
        to_plot.areas_info(iarea).name);
    saveas(fig,outfile);
    close(fig);
    
    % Select broadban frequencies
    f_index = to_plot.bands_info(6).f_limits_index;
    f_of_interest = to_plot.bands_info(6).f_original;
    f_of_interest = f_of_interest(f_index);
    
    % We are plotting a figure for each tester vs sEEGnal
    for itester = 1 : size(pow_human_dataset,4)
        
        fig = figure('WindowState','maximized');
        hold on
        
        current_pow_human_dataset = pow_human_dataset(channels_index,f_index,:,itester);
        current_pow_sEEGnal_dataset = pow_sEEGnal_dataset(channels_index,f_index,:);
        
        % Average across electrodes
        human_average = squeeze(nanmean(nanmean(current_pow_human_dataset,1),3));
        sEEGnal_average = squeeze(nanmean(nanmean(current_pow_sEEGnal_dataset,1),3));
        
        % Plot  
        plot(f_of_interest,squeeze(nanmean(current_pow_human_dataset,1)),...
            'Color',[0.6980    0.8902    0.8902],'LineWidth',0.5)
        plot(f_of_interest,squeeze(nanmean(current_pow_sEEGnal_dataset,1)),...
            'Color',[1.0000    0.7725    0.6902],'LineWidth',0.5)
        plot(f_of_interest,human_average,'Color',[0 0.5 0.5],'LineWidth',3)
        plot(f_of_interest,sEEGnal_average,'Color',[1 0.498 0.314],'LineWidth',3)
        
        % Enhance the plot
        lines = findobj(gca, 'Type', 'Line');
        legend(lines(1:2), {'sEEGnal',config.testers{itester}});
        title(sprintf('Power in BroadBand and %s area',  to_plot.areas_info(iarea).name),...
            'Interpreter','none')
        
        % Save the figure
        outfile = sprintf('%s/%s_%s_pow_spectrum.svg',config.path.figures,...
            config.testers{itester},to_plot.areas_info(iarea).name);
        saveas(fig,outfile);
        close(fig);
        
    end
    
    
    
    % We are plotting a figure with the average across channels and then
    % across recordings for all testers vs sEEGnal
    % First, plot the testers
    fig = figure('WindowState','maximized');
    hold on
    colors_mean = [0 0.4470 0.7410;
        0.3010 0.3 1;
        0.4940 0.1840 0.5560;
        0.6350 0.0780 0.1840];
    for itester = 1 : size(pow_human_dataset,4)
        
        current_pow_human_dataset = pow_human_dataset(channels_index,f_index,:,itester);
        human_average = squeeze(nanmean(nanmean(current_pow_human_dataset,1),3));
        
        % Plot
        plot(f_of_interest,human_average,...
            'Color',colors_mean(itester,:),'LineWidth',3)
        
    end
    
    % Now plot sEEGnal
    sEEGnal_average = squeeze(nanmean(nanmean(current_pow_sEEGnal_dataset,1),3));    
    plot(f_of_interest,sEEGnal_average,'Color',[1 0.498 0.314],'LineWidth',3)
    
    % Enhance the plot
    legend(cat(2,config.testers,'sEEGnal'));
    title(sprintf('Power in BroadBand and %s area',  to_plot.areas_info(iarea).name),...
        'Interpreter','none')
    
    % Save the figure
    outfile = sprintf('%s/all_%s_pow_spectrum.svg',config.path.figures,...
        to_plot.areas_info(iarea).name);
    saveas(fig,outfile);
    close(fig);
    
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
% pos_labels = [pos_elec(:,1)-0.015, pos_elec(:,2)+0.03];
% text(pos_labels(:,1),pos_labels(:,2),config.complete_channel_labels)

end