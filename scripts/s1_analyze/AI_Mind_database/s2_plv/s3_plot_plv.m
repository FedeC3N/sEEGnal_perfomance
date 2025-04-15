%{

Different plots to investigate the results.


@author: Fede

%}

clear
clc
close all
restoredefaultpath

% Add paths
addpath('../shared/')

% Paths
config.path.clean_data = '../../../../databases/AI_Mind_database/derivatives';
config.path.results = '../../../../results/AI_Mind_database/plv';
config.path.figures = '../../../../docs/manuscript/figures/AI_Mind_database/plv_results';

if ~exist(config.path.figures), mkdir(config.path.figures),end

% Load the results
load(sprintf('%s/plv_results.mat',config.path.results));

% Get the different testers (except sEEGnal)
testers = dir(sprintf('%s/*',config.path.clean_data));
testers = testers(3:end-1);
testers = {testers.name};
config.testers = testers;

% Define the frequency bands
bands_info = struct('name',{'delta', 'theta','alpha','beta','gamma'},...
    'f_limits',{[ 2 4] , [4 8] , [8 12] , [12 20] , [20 45] });

% Areas
% Areas
areas_info = struct('name',{'frontal','temporal_l','temporal_r','parietal_l',...
    'parietal_r','occipital','sorted'},'channel',[]);
areas_info(1).channel = {'Fp1', 'Fp2', 'F7', 'F3', 'Fz', 'F4', 'F8', 'AF7', 'AF3', 'AF4', 'AF8', 'F5', 'F1', 'F2', 'F6'};
areas_info(2).channel = {'T7', 'FT7', 'TP7'};
areas_info(3).channel = {'T8', 'FT8', 'TP8'};
areas_info(4).channel = {'CP5', 'CP1', 'P7', 'P3', 'TP7', 'CP3', 'P5', 'P1'};
areas_info(5).channel = {'CP2', 'CP6', 'P4', 'P8', 'CP4', 'TP8', 'P2', 'P6'};
areas_info(6).channel = {'PO9', 'O1', 'Oz', 'O2', 'PO10','PO7', 'PO3', 'POz', 'PO4', 'PO8'};
areas_info(7).channel = {'Fp1', 'Fp2', 'F7', 'F3', 'Fz', 'F4', 'F8', 'FC5', 'FC1', 'FC2', 'FC6', 'AFz','AF7', 'AF3', 'AF4', 'AF8', 'F5', 'F1', 'F2', 'F6', 'FT7', 'FC3', 'FC4', 'FT8',...
    'T7','T8','TP7','TP8',...
    'C3', 'Cz', 'C4', 'CP5', 'CP1', 'CP2', 'CP6','C5', 'C1', 'C2', 'C6','CP3', 'CPz', 'CP4',...
    'P7', 'P3', 'Pz', 'P4', 'P8', 'PO9','PO10','P5', 'P1', 'P2', 'P6', 'PO7', 'PO3', 'POz', 'PO4', 'PO8',...
    'O1', 'Oz', 'O2'};

% Define the frequency bands
config.bands = {'delta', 'theta','alpha','beta','gamma'};

% Define measures
config.measures = {'NRMSE', 'rho'};

% Channels
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

% Read the power spectrum
[plv_human_dataset] = read_plv_dataset(config,testers);
[plv_sEEGnal_dataset] = read_plv_dataset(config,{'sEEGnal'});

%%%%%%%%%%%%%%%
% PLOTS
%%%%%%%%%%%%%%%

% Topoplots of the stats
plot_stats_in_head(config,stats)

% Plot differences in channels as measures by NRMSE
plot_NRMSE_violinplots(config,stats,bands_info)

% Plot differences in channels as measures by corr
plot_corr(config,bands_info,plv_human_dataset,plv_sEEGnal_dataset)


% Aux functions
function plv_dataset = read_plv_dataset(config,dataset_name)

for itester = 1 : numel(dataset_name)
    
    % Load the datset
    dataset_path = sprintf('%s/%s/%s_dataset.mat',config.path.clean_data,...
        dataset_name{itester},dataset_name{itester});
    dummy = load(dataset_path);
    
    for icurrent = 1 : numel(dummy.dataset)
        
        % Load pow
        current_dataset = dummy.dataset(icurrent);
        plv_file = sprintf('../../../../%s/%s.mat',current_dataset.plv.path,...
            current_dataset.plv.file);
        
        if ~exist(plv_file)
            fake_pow = nan(size(current_pow));
            pow_dataset_norm(:,:,icurrent,itester) = fake_pow;
            continue
        end
        
        plv = load(plv_file);
        
        % Create the PLV_all struct and the channels info
        if icurrent == 1
            n_sensors = numel(plv.channels_included_index);
            n_bands = numel(plv.bands_info);
            n_subjects = numel(dummy.dataset);
            n_testers = numel(dataset_name);
            plv_dataset = nan(n_sensors,n_sensors,...
                n_bands,n_subjects,n_testers);

        end
        
        % Save each PLV
        for iband = 1 : numel(plv.bands_info)
            
            plv_dataset(:,:,iband,icurrent,itester) = plv.(plv.bands_info(iband).name).plv;
            
        end
            
    end
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
        ax1 = subplot(1,2,imeasure);
        [pos_elec, size_elec] = draw_head(config);
        
        % Scatter the values of interest
        current_stats = nanmean(stats.(current_measure)(:,:,iband,:),4);
        color_elec = nanmean(current_stats,2);
        scatter(pos_elec(:,1),pos_elec(:,2),size_elec,color_elec,'filled')
        c = colorbar;
        c.Location = 'south';
           
        % Set the colormap
        customCMap = [73 144 209; 133 58 123; 209 70 2]/255; 
        colormap(customCMap)
        
        % Details
        title(sprintf('Band %s - Measure %s', current_band, current_measure))
           
    end
    
    % Save the figure
    outfile = sprintf('%s/%s_head_measures.svg',config.path.figures,...
        current_band);
    saveas(fig,outfile);
    close(fig);
    
    
end

end


function plot_NRMSE_violinplots(config,stats,bands_info)

fig = figure('WindowState', 'maximized');
hold on

% For each band
for iband = 1 : numel(bands_info)
    
    current_band = bands_info(iband).name;
    current_NRMSE = stats.NRMSE(:,:,iband,:);
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
legend(bands_info.name)

% Save the figure
outfile = sprintf('%s/NRMSE_violinplot.svg',config.path.figures);
saveas(fig,outfile);
close(fig);

end


function plot_corr(config,bands_info,plv_human_dataset,plv_sEEGnal_dataset)

colors = [0    0.4470    0.7410;...
    0.8500    0.3250    0.0980;...
    0.9290    0.6940    0.1250;...
    0.4940    0.1840    0.5560;...
    0.4660    0.6740    0.1880];

% For each band
fig = figure('WindowState', 'maximized');
hold on
for iband = 1 : numel(bands_info)
    
 
    % Get the pow in the current frequencies
    current_band = bands_info(iband).name;
    current_human = squeeze(plv_human_dataset(:,:,iband,:,:));
    current_human = nanmean(current_human,4);
    % Create the upper matrix mask
    triu_mask = triu(ones(size(current_human(:,:,1))));
    triu_mask = repmat(triu_mask,1,1,size(current_human,3));
    triu_mask = logical(triu_mask);
    current_human = current_human(triu_mask);
    current_sEEGnal = squeeze(plv_sEEGnal_dataset(:,:,iband,:));
    current_sEEGnal = current_sEEGnal(triu_mask);
    
    % Remove the nans
    nans_mask = isnan(current_human) | isnan(current_sEEGnal);
    current_human = current_human(~nans_mask);
    current_sEEGnal = current_sEEGnal(~nans_mask);
    
    % Plot
    subplot(2,3,iband)
    s = scatter(current_human,current_sEEGnal,'filled');
    s.SizeData = 10;
    s.MarkerEdgeColor = colors(iband,:);
    s.MarkerFaceColor = colors(iband,:);
    s.MarkerFaceAlpha = 0.2;
    s.MarkerEdgeAlpha = s.MarkerFaceAlpha;
    
    % Set the limits of the axis to be square
    lims = axis;
    xlim([min(lims) max(lims)])
    ylim([min(lims) max(lims)])
    axis square
    
    % lsline
    ls = lsline;
    ls.Color = [255 136 136]/255;
    ls.LineWidth = 2;
    
    % Plot the identity line
    l = line([min(lims) max(lims)],[min(lims) max(lims)],'Color','black','LineStyle','--');
    
    % Info
    title(sprintf('%s band correlation',current_band))
    xlabel('Human Expert')
    ylabel('sEEGnal')

end

    
% Save the figure
outfile = sprintf('%s/plv_corr.svg',config.path.figures);
saveas(fig,outfile);
close(fig);

end

