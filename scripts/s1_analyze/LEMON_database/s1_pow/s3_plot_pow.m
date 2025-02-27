clear
clc
close all
restoredefaultpath

% Paths
config.path.clean_data = '../../../../databases/LEMON_database/derivatives';
config.path.results = '../../../../results/LEMON_database/pow';
config.path.figures = '../../../../docs/manuscript/figures/LEMON_database/pow_results';

% Add paths
addpath('../shared/')

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
plot_pow_spectrum_norm_one_subject(config,pow_lemon_dataset,pow_sEEGnal_dataset)

% Violinplots
plot_NRMSE_violinplots(config,stats,bands_info)

% corr
plot_corr(config,bands_info,pow_lemon_dataset,pow_sEEGnal_dataset)

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

for iband = 1 : numel(config.bands_info)
    
    current_band = config.bands_info(iband).name;
    
    % Define figure
    fig = figure('WindowState', 'maximized');
    for imeasure = 1 : numel(config.measures)
        
        current_measure = config.measures{imeasure};
        
        % Scatter the head
        ax1 = subplot(1,2,imeasure);
        [pos_elec, size_elec] = draw_head(config);
        current_stats = stats.(current_band).(current_measure);
        color_elec = nanmean(current_stats,2);
        
        % Scatter the values of interest
        scatter(pos_elec(:,1),pos_elec(:,2),size_elec,color_elec,'filled')
        c = colorbar;
        c.Location = 'south';
        
        % Set the colormap
        customCMap = [73 144 209; 133 58 123; 209 70 2]/255; 
        colormap(customCMap)
        
        % Details
        title(sprintf('%s band - Measure %s',current_band, current_measure))
        
        
    end
    
    % Save the figure
    outfile = sprintf('%s/%s_whole_head_head_measures.svg',...
        config.path.figures, current_band);
    saveas(fig,outfile);
    close(fig);
end

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


function plot_NRMSE_violinplots(config,stats,bands_info)

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
legend(bands_info.name)

% Save the figure
outfile = sprintf('%s/pow_NRMSE_violinplot.svg',config.path.figures);
saveas(fig,outfile);
close(fig);


end


function plot_corr(config,bands_info,pow_lemon_dataset,pow_sEEGnal_dataset)

% Remove the broadband
bands_info = bands_info(1:end-1);

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
    current_f = bands_info(iband).f_limits_index;
    current_human = pow_lemon_dataset(:,current_f,:);
    current_human = current_human(:);
    current_sEEGnal = pow_sEEGnal_dataset(:,current_f,:);
    current_sEEGnal = current_sEEGnal(:);
    
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
outfile = sprintf('%s/pow_corr.svg',config.path.figures);
saveas(fig,outfile);
close(fig);


end

