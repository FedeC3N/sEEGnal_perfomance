clear
clc
close all
restoredefaultpath

% Paths
config.path.clean_data = '../../../../databases/LEMON_database/derivatives';
config.path.results = '../../../../results/LEMON_database/pow';
config.path.figures = '../../../../docs/manuscript/figures/LEMON_database/pow_results';

% Load the results
load(sprintf('%s/pow_results.mat',config.path.results));

% Get the different testers
testers = {'lemon','sEEGnal'};

% Channels
config.complete_channel_labels = {'Fp1', 'Fp2', 'F7', 'F3', 'Fz', 'F4', 'F8', 'FC5', 'FC1', 'FC2', 'FC6', 'T7', 'C3', 'Cz', 'C4', 'T8', 'CP5', 'CP1', 'CP2', 'CP6', 'AFz', 'P7', 'P3', 'Pz', 'P4', 'P8', 'PO9', 'O1', 'Oz', 'O2', 'PO10', 'AF7', 'AF3', 'AF4', 'AF8', 'F5', 'F1', 'F2', 'F6', 'FT7', 'FC3', 'FC4', 'FT8', 'C5', 'C1', 'C2', 'C6', 'TP7', 'CP3', 'CPz', 'CP4', 'TP8', 'P5', 'P1', 'P2', 'P6', 'PO7', 'PO3', 'POz', 'PO4', 'PO8'};

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
config.bands_info = bands_info;

% Define measures
% config.measures = {'NRMSE', 'rho', 'tstat'};
config.measures = {'NRMSE', 'rho'};

% Read the power spectrum
pow_lemon_dataset = read_pow_dataset(config,{'lemon'});
pow_sEEGnal_dataset = read_pow_dataset(config,{'sEEGnal'});

%%%%%%%%%%%%%%%
% PLOTS
%%%%%%%%%%%%%%%

% corr and NRMSE in head
plot_stats_in_head(config,stats)

% Normalized pow spectrum for all areas in broadband
% Prepare the pow_spectrum
plot_pow_spectrum_norm_one_subject(config,pow_lemon_dataset,pow_sEEGnal_dataset);

% Violinplots
plot_NRMSE_channels(config,stats,bands_info)

% Functions
function pow_dataset_norm = read_pow_dataset(config,dataset_name)

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


function plot_stats_in_head(config,stats)

% Define figure
fig = figure('WindowState', 'maximized');
for imeasure = 1 : numel(config.measures)
    
    current_measure = config.measures{imeasure};
    
    % Scatter the head
    ax1 = subplot(1,2,imeasure);
    [pos_elec, size_elec] = draw_head(config);
    current_stats = stats.broadband.(current_measure);
    color_elec = nanmean(current_stats,2);
    
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
outfile = sprintf('%s/whole_head_head_measures.svg',config.path.figures);
saveas(fig,outfile);
close(fig);
    

end


function plot_pow_spectrum_norm_one_subject(config,pow_human_dataset,pow_sEEGnal_dataset)


% Select one subject to plot as an example
% I checked all the recordings and 27 is the best picture
for random_index = 27
    
    % Estimate mean error
    fig = figure('WindowState','maximized');
    hold on
    f_of_interest = config.bands_info(6).f_original;
    
    % Estimate the powers
    current_pow_human_dataset = pow_human_dataset(:,:,random_index);
    human_average = nanmean(nanmean(current_pow_human_dataset,3),1);
    current_pow_sEEGnal_dataset = pow_sEEGnal_dataset(:,:,random_index);
    sEEGnal_average = nanmean(current_pow_sEEGnal_dataset,1);
    
    % Plot channels first
    plot(f_of_interest,nanmean(current_pow_human_dataset,3),...
        'Color',[0.6980    0.8902    0.8902],'LineWidth',0.5)
    plot(f_of_interest,nanmean(current_pow_sEEGnal_dataset,4),...
        'Color',[1.0000    0.7725    0.6902],'LineWidth',0.5)
    
    % Plot mean
    plot(f_of_interest,human_average,...
        'Color',[0 0.5 0.5],'LineWidth',3)
    plot(f_of_interest,sEEGnal_average,...
        'Color',[1 0.498 0.314],'LineWidth',3)
    
    % Enhance the plot
    lines = findobj(gca, 'Type', 'Line');
    legend(lines(1:2),{'sEEGnal','Human experts'});
    title('Power in BroadBand in whole head')
    
    % Save the figure
    outfile = sprintf('%s/%i_pow_spectrum.svg',config.path.figures,...
        random_index);
    saveas(fig,outfile);
    close(fig);
    
end

end


function plot_NRMSE_channels(config,stats,bands_info)

fig = figure('WindowState', 'maximized');
hold on

% Remove the broadband
bands_info = bands_info(1:end-1);

% For each band
for iband = 1 : numel(bands_info)
    
    current_band = bands_info(iband).name;
    current_NRMSE = stats.(current_band).NRMSE;
    current_NRMSE = current_NRMSE(:);
    
    % X axis for plot
    x_vector = iband * ones(numel(current_NRMSE),1);
    
    % Plot
    sw = swarmchart(x_vector,current_NRMSE,'filled','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5);
    bx = boxchart(x_vector,current_NRMSE,'BoxFaceColor',sw.CData ./ 1.2,'WhiskerLineColor',sw.CData ./ 1.2,...
        'MarkerStyle','none','BoxWidth',sw.XJitterWidth);
    
    
end

xlim([0 numel(bands_info) + 1])
xticks(1:numel(bands_info))
xticklabels({bands_info.name})
title('NRMSE for each channel')
ylabel('NRMSE','Interpreter','none')
set(gca,'TickLabelInterpreter','none')

% Save the figure
outfile = sprintf('%s/pow_NRMSE_violinplot.svg',config.path.figures);
saveas(fig,outfile);
close(fig);


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