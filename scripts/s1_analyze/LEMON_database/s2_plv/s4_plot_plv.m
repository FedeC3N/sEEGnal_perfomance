clear
clc
close all
restoredefaultpath

% Add paths
addpath('../shared/')

% Paths
config.path.clean_data = '../../../../databases/LEMON_database/derivatives';
config.path.results = '../../../../results/plv';

% Load the results
load(sprintf('%s/plv_results.mat',config.path.results));

% Get the different testers
testers = {'lemon','sEEGnal'};

% Read the power spectrum
[plv_lemon_dataset,~,channels_lemon_included] = read_plv_dataset(config,'lemon');
[plv_sEEGnal_dataset,f,channels_sEEGnal_included] = read_plv_dataset(config,'sEEGnal');

% Define the frequency bands
bands_info = struct('name',{'delta', 'theta','alpha','beta','gamma'},...
    'f_limits',{[ 2 4] , [4 8] , [8 12] , [12 20] , [20 45] });

% Areas
areas_info = struct('name',{'frontal','temporal_l','temporal_r','parietal_l',...
    'parietal_r','occipital'},'channel',[]);
areas_info(1).channel = {'Fp1';'AF7';'AF3'; 'Fpz';'Fp2';'AF8';'AF4';'AFz'};
areas_info(2).channel = {'FT7';'FT9';'T7';'TP7';'TP9'};
areas_info(3).channel = {'FT8';'FT10';'T8';'TP8';'TP10'};
areas_info(4).channel = {'CP5';'CP3';'CP1';'P1';'P3';'P5';'P7';'P9'};
areas_info(5).channel = {'CP6';'CP4';'CP2';'P2';'P4';'P6';'P8';'P10'};
areas_info(6).channel = {'O1';'Iz';'Oz';'O9';'O2';'O10'};

% Define the frequency bands
config.bands = {'delta', 'theta','alpha','beta','gamma'};

% Define measures
config.measures = {'NMSE', 'rho', 'tstat'};

% Channels
config.complete_channel_labels = {'Fp1', 'Fp2', 'F7', 'F3', 'Fz', 'F4', 'F8', 'FC5', 'FC1', 'FC2', 'FC6', 'T7', 'C3', 'Cz', 'C4', 'T8', 'CP5', 'CP1', 'CP2', 'CP6', 'AFz', 'P7', 'P3', 'Pz', 'P4', 'P8', 'PO9', 'O1', 'Oz', 'O2', 'PO10', 'AF7', 'AF3', 'AF4', 'AF8', 'F5', 'F1', 'F2', 'F6', 'FT7', 'FC3', 'FC4', 'FT8', 'C5', 'C1', 'C2', 'C6', 'TP7', 'CP3', 'CPz', 'CP4', 'TP8', 'P5', 'P1', 'P2', 'P6', 'PO7', 'PO3', 'POz', 'PO4', 'PO8'};

%%%%%%%%%%%%%%%
% PLOTS
%%%%%%%%%%%%%%%

% Topoplots of the stats
plot_stats_in_head(config,stats)

% Plot difference in channels
plot_diff_channels(plv_lemon_dataset,plv_sEEGnal_dataset,bands_info)

% Plot correlation of each pair of channels
plot_corr_channels(plv_lemon_dataset,plv_sEEGnal_dataset,bands_info, channels_sEEGnal_included)

% Plot correlation of PLV matrix for each subject
plot_corr_plv_subjects(plv_lemon_dataset,plv_sEEGnal_dataset,bands_info, channels_lemon_included,channels_sEEGnal_included)

% Plot PLV correlation by areas
plot_plv_corr_areas(config,plv_lemon_dataset,plv_sEEGnal_dataset,bands_info,areas_info)

% Aux functions
function [plv_dataset,bands_info,channels] = read_plv_dataset(config, dataset_name)

% Load the datset
dataset_path = sprintf('%s/%s/%s_dataset.mat',config.path.clean_data,...
    dataset_name,dataset_name);
dummy = load(dataset_path);

for icurrent = 1 : numel(dummy.dataset)
    
    % Load plv
    plv = load(sprintf('../../../../%s/%s',dummy.dataset(icurrent).plv.path,...
        dummy.dataset(icurrent).plv.file));
    
    % Save bands_info
    bands_info = plv.bands_info;
    
    % Create the PLV_all struct and the channels info
    if icurrent == 1
        n_sensors = numel(plv.channels_included_index);
        n_bands = numel(bands_info);
        n_subjects = numel(dummy.dataset);
        plv_dataset = nan(n_sensors,n_sensors,...
            n_bands,n_subjects);
        channels = struct('channels_included',[],...
            'channels_included_index',[]);
    end
    
    % Save each PLV
    for iband = 1 : numel(bands_info)
        
        plv_dataset(:,:,iband,icurrent) = plv.(bands_info(iband).name).plv;
        
    end
    
    % Save the channel information
    channels(icurrent).channels_included = plv.channels_included;
    channels(icurrent).channels_included_index = plv.channels_included_index;
    
end

end


function plot_stats_in_head(config,stats)

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
        if min(color_elec) < 0
            caxis([min(color_elec), max(color_elec)]);
        else
            caxis([0, max(color_elec)]);
        end
        
        % Set the colormap
        gradient_1_neg = linspace(1,1,200);
        gradient_2_neg = linspace(1,0.4,200);
        gradient_3_neg = linspace(1,0.4,200);
        gradient_1_pos = linspace(0.4,1,200);
        gradient_2_pos = linspace(0.4,1,200);
        gradient_3_pos = linspace(1,1,200);
        clrmap = [[gradient_1_pos gradient_1_neg]',...
            [gradient_2_pos gradient_2_neg]',...
            [gradient_3_pos gradient_3_neg]'];
        if min(color_elec) > 0
            clrmap = clrmap(200:end,:);
        end
        colormap(ax1,clrmap)
        
        
        
        
        % Details
        title(sprintf('Band %s - Measure %s', current_band, current_measure))
        
        
    end
    
    
end

end


function plot_diff_channels(plv_lemon_dataset,plv_sEEGnal_dataset,bands_info)

figure('WindowState', 'maximized');
hold on

% For each band
for iband = 1 : numel(bands_info)
    
    % Select the band data
    current_plv_global_lemon = squeeze(plv_lemon_dataset(:,:,iband,:));
    current_plv_global_sEEGnal = squeeze(plv_sEEGnal_dataset(:,:,iband,:));
    
    % Get the difference of each channel and average all the differences
    plv_diff = current_plv_global_lemon - current_plv_global_sEEGnal;
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
title('Average PLV difference (in %) for each channel')
ylabel('% Variation','Interpreter','none')
set(gca,'TickLabelInterpreter','none')



end


function plot_corr_channels(plv_lemon_dataset,plv_sEEGnal_dataset,bands_info, channels)

figure('WindowState', 'maximized');
hold on

% For each band
for iband = 1 : numel(bands_info)
    
    % Select the band data
    current_plv_global_lemon = squeeze(plv_lemon_dataset(:,:,iband,:));
    current_plv_global_sEEGnal = squeeze(plv_sEEGnal_dataset(:,:,iband,:));
    plv_vector_lemon = reshape(current_plv_global_lemon,[],size(current_plv_global_lemon,3));
    plv_vector_sEEGnal = reshape(current_plv_global_sEEGnal,[],size(current_plv_global_sEEGnal,3));
    
    % Get the upper matrix
    dummy = ones(size(current_plv_global_sEEGnal(:,:,1)));
    upper_index = abs(triu(dummy,1)) > 0 ;
    upper_index = upper_index(:);
    plv_vector_lemon =  plv_vector_lemon(upper_index,:);
    plv_vector_sEEGnal = plv_vector_sEEGnal(upper_index,:);
    
    % Estimate correlation for each channel
    plv_vector_lemon = plv_vector_lemon';
    plv_vector_sEEGnal = plv_vector_sEEGnal';
    rho = nan(1,size(plv_vector_sEEGnal,2));
    x_vector = iband * ones(1,size(plv_vector_sEEGnal,2));
    for ichannel = 1 : size(plv_vector_sEEGnal,2)
        
        % Remove nans
        current_channel_lemon = plv_vector_lemon(:,ichannel);
        current_channel_sEEGnal = plv_vector_sEEGnal(:,ichannel);
        nan_index = isnan(current_channel_lemon) | isnan(current_channel_sEEGnal);
        current_channel_lemon = current_channel_lemon(~nan_index);
        current_channel_sEEGnal = current_channel_sEEGnal(~nan_index);
        
        
        rho(ichannel) = corr(current_channel_lemon,current_channel_sEEGnal);
        
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


function plot_corr_plv_subjects(plv_lemon_dataset,plv_sEEGnal_dataset,bands_info, channels_lemon,channels_sEEGnal)

figure('WindowState', 'maximized');
hold on

% For each band
for iband = 1 : numel(bands_info)
    
    % Create the rho empty and the plot axis
    rho = nan(1,size(plv_lemon_dataset,4));
    x_vector = iband * ones(1,size(rho,2));
    for irho = 1 :numel(rho)
        
        % Get the current plv matrix
        current_plv_lemon = plv_lemon_dataset(:,:,iband,irho);
        current_plv_sEEGnal = plv_sEEGnal_dataset(:,:,iband,irho);
        
        % Create an index of valid channels
        valid = channels_lemon(irho).channels_included_index & channels_sEEGnal(irho).channels_included_index;
        
        % Remove badchannels
        current_plv_lemon = current_plv_lemon(valid,valid);
        current_plv_sEEGnal = current_plv_sEEGnal(valid,valid);
        
        % Get the upper matrix in vector format
        dummy = ones(size(current_plv_sEEGnal));
        upper_index = abs(triu(dummy,1)) > 0 ;
        current_plv_lemon = current_plv_lemon(upper_index);
        current_plv_sEEGnal = current_plv_sEEGnal(upper_index);
        
        % Correlation
        rho(irho) = corr(current_plv_lemon,current_plv_sEEGnal);
        
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


function plot_plv_corr_areas(config,plv_lemon_dataset,plv_sEEGnal_dataset,bands_info,areas_info)

% Define figure
fig = figure('WindowState', 'maximized');

% For each band
for iband = 1 : numel (bands_info)
    
    % Select the matrix
    current_plv_lemon = squeeze(plv_lemon_dataset(:,:,iband,:));
    current_plv_sEEGnal = squeeze(plv_sEEGnal_dataset(:,:,iband,:));
    
    % For each area
    corr_matrix = nan(numel(areas_info));
    NRMSE_matrix = nan(numel(areas_info));
    for iarea_1 = 1 : numel(areas_info)
        
        % Select the matrix
        channels_mask = ismember(config.complete_channel_labels,areas_info(iarea_1).channel);
        channel_index = find(channels_mask);
        current_plv_lemon_area_1 = current_plv_lemon(:,channels_mask,:);
        current_plv_sEEGnal_area_1 = current_plv_sEEGnal(:,channels_mask,:);
        
        for iarea_2 = 1 : numel(areas_info)
            
            % Select the matrix
            channels_mask = ismember(config.complete_channel_labels,areas_info(iarea_2).channel);
            channel_index = find(channels_mask);
            current_plv_lemon_area_2 = current_plv_lemon_area_1(channels_mask,:,:);
            current_plv_sEEGnal_area_2 = current_plv_sEEGnal_area_1(channels_mask,:,:);
            
            % Create a corr and NRMSE matrix
            dummy_corr_matrix = nan(size(current_plv_lemon_area_2,2),size(current_plv_lemon_area_2,3));
            dummy_NRMSE = nan(size(current_plv_lemon_area_2,2),size(current_plv_lemon_area_2,3));
            
            % For each channel
            for ichannel = 1 : size(current_plv_lemon_area_2,2)
                
                % Get the current channel and the rest of the matrix
                current_channel_lemon = squeeze(current_plv_lemon_area_2(:,ichannel,:));
                current_channel_sEEGnal = squeeze(current_plv_sEEGnal_area_2(:,ichannel,:));
                
                % Estimate NRMSE and corr
                for isubject = 1 : size(dummy_corr_matrix,2)
                    
                    dummy_corr_matrix(ichannel,isubject) = corr(current_channel_lemon(:,isubject),current_channel_sEEGnal(:,isubject),'Rows','complete');
                    
                end
                
            end
            
            % Average across subjects
            dummy_corr_matrix = nanmean(dummy_corr_matrix,2);
            
            % Average across channels and save in the matrix
            corr_matrix(iarea_1,iarea_2) = nanmean(dummy_corr_matrix);
            
        end
        
        
    end
    
    % Plot
    ax1 = subplot(1,numel(bands_info),iband);
    imagesc(corr_matrix)
    axis square
    
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
    if min(corr_matrix(:)) > 0
        clrmap = clrmap(2000:end,:);
    end
    colormap(ax1,clrmap)
    colorbar
    
end


end
