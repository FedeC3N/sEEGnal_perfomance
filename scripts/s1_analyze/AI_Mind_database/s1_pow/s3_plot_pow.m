%{

Different plots to investigate the results.


@author: Fede

%}

clear
clc
close all
restoredefaultpath

% Paths
config.path.clean_data = '../../../../databases/AI_Mind_database/derivatives';
config.path.results = '../../../../results/AI_Mind_database/pow';
config.path.figures = '../../../../docs/manuscript/figures/AI_Mind_database/pow_results';
if ~exist(config.path.figures), mkdir(config.path.figures), end

% Add paths
addpath('../shared/')

% Load the results
load(sprintf('%s/pow_results.mat',config.path.results));

% Get the different testers (except sEEGnal)
testers = dir(sprintf('%s/*',config.path.clean_data));
testers = testers(3:end-1);
testers = {testers.name};
config.testers = testers;

% Channels
config.complete_channel_labels = {'Fp1', 'Fpz', 'Fp2', 'F7', 'F3', 'Fz', 'F4', 'F8', 'FC5', 'FC1', 'FC2', 'FC6',...
    'T7', 'C3', 'Cz', 'C4', 'T8', 'CP5', 'CP1', 'CP2', 'CP6', 'P7', 'P3',...
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
config.bands_info = bands_info;

% Define measures
% config.measures = {'NRMSE', 'rho', 'tstat'};
config.measures = {'NRMSE'};

% Read the power spectrum
[pow_human_dataset] = read_pow_dataset(config,testers);
[pow_sEEGnal_dataset] = read_pow_dataset(config,{'sEEGnal'});

%%%%%%%%%%%%%%%
% PLOTS
%%%%%%%%%%%%%%%

% corr and NRMSE in head
% plot_stats_in_head(config,stats)

% Normalized pow spectrum for all areas in broadband
% plot_pow_spectrum_norm_one_subject(config,pow_human_dataset,pow_sEEGnal_dataset);

% Violinplots
% plot_violinplots(config,stats,bands_info)

% corr
plot_sEEGnal_vs_human(config,bands_info,pow_human_dataset,pow_sEEGnal_dataset)

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



function plot_stats_in_head(config,stats)

for imeasure = 1 : numel(config.measures)

    current_measure = config.measures{imeasure};

    % Colorbar limits
    dummy = nanmean(stats.(current_measure),2);
    high_limit = max(dummy(:));

    % Define figure
    fig = figure('WindowState', 'maximized');

    for iband = 1 : numel(config.bands_info)

        current_band = config.bands_info(iband).name;

        % Scatter the head
        ax1 = subplot(2,3,iband);
        [pos_elec, size_elec] = draw_head(config);
        current_stats = nanmean(stats.(current_measure)(:,:,iband,:),4);
        color_elec = nanmean(current_stats,2);

        % Scatter the values of interest
        scatter(pos_elec(:,1),pos_elec(:,2),size_elec,color_elec,'filled')
        c = colorbar;
        c.Location = 'south';

        % Set the colormap
        switch current_measure
            case 'NRMSE'
                caxis([0, high_limit])
                customCMap = 'jet';
            case 'rho'
                caxis([-1 1]);
                customCMap = 'jet';
        end
        colormap(ax1,customCMap)

        % Details
        title(sprintf('%s band',current_band))


    end

    % Details
    sgtitle(sprintf('Measure %s', current_measure))

    % Save the figure
    outfile = sprintf('%s/%s_whole_head.svg',...
        config.path.figures,current_measure);
    saveas(fig,outfile);
    outfile = sprintf('%s/%s_whole_head.png',...
        config.path.figures,current_measure);
    saveas(fig,outfile);
    close(fig);

end

end



function plot_pow_spectrum_norm_one_subject(config,pow_human_dataset,pow_sEEGnal_dataset)


% Select one subject to plot as an example
% I checked all the recordings and 19 is the best picture
for random_index = 19

    % Estimate mean error
    fig = figure('WindowState','maximized');
    hold on
    f_of_interest = config.bands_info(6).f_original;

    % Estimate the powers
    current_pow_human_dataset = squeeze(pow_human_dataset(:,:,random_index,:));
    human_average = nanmean(nanmean(current_pow_human_dataset,3),1);
    current_pow_sEEGnal_dataset = pow_sEEGnal_dataset(:,:,random_index);
    sEEGnal_average = nanmean(current_pow_sEEGnal_dataset,1);

    % Plot channels first
    plot_sEEGnal_vs_human(f_of_interest,nanmean(current_pow_human_dataset,3),...
        'Color',[0.6980    0.8902    0.8902],'LineWidth',0.5)
    plot_sEEGnal_vs_human(f_of_interest,nanmean(current_pow_sEEGnal_dataset,4),...
        'Color',[1.0000    0.7725    0.6902],'LineWidth',0.5)

    % Plot mean
    plot_sEEGnal_vs_human(f_of_interest,human_average,...
        'Color',[0 0.5 0.5],'LineWidth',3)
    plot_sEEGnal_vs_human(f_of_interest,sEEGnal_average,...
        'Color',[1 0.498 0.314],'LineWidth',3)

    % Enhance the plot
    lines = findobj(gca, 'Type', 'Line');
    legend(lines(1:2),{'sEEGnal','Human experts'});
    title('Power in BroadBand in whole head')

    % Save the figure
    outfile = sprintf('%s/%i_pow_spectrum.svg',config.path.figures,...
        random_index);
    saveas(fig,outfile);
    outfile = sprintf('%s/%i_pow_spectrum.png',config.path.figures,...
        random_index);
    saveas(fig,outfile);
    close(fig);

end

end



function plot_violinplots(config,stats,bands_info)

% Remove the broadband
bands_info = bands_info(1:end-1);

for imeasure = 1 : numel(config.measures)

    current_measure = config.measures{imeasure};

    fig = figure('WindowState', 'maximized');
    hold on

    % For each band
    for iband = 1 : numel(bands_info)

        current_band = bands_info(iband).name;
        current_stat = stats.(current_measure);
        current_stat = current_stat(:,:,iband,:);
        current_stat = current_stat(:);

        % X axis for plot
        x_vector = iband * ones(numel(current_stat),1);

        % Plot
        sw = swarmchart(x_vector,current_stat,'filled','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5);
        bx = boxchart(x_vector,current_stat,'BoxFaceColor',sw.CData ./ 1.2,'WhiskerLineColor',sw.CData ./ 1.2,...
            'MarkerStyle','none','BoxWidth',sw.XJitterWidth);


    end

    xlim([0 numel(bands_info) + 1])
    xticks(1:numel(bands_info))
    xticklabels({bands_info.name})
    title(sprintf('%s for each channel',current_measure))
    ylabel(current_measure,'Interpreter','none')
    set(gca,'TickLabelInterpreter','none')
    legend(bands_info.name)

    % Save the figure
    outfile = sprintf('%s/pow_%s_violinplot.svg',config.path.figures,current_measure);
    saveas(fig,outfile);
    outfile = sprintf('%s/pow_%s_violinplot.png',config.path.figures,current_measure);
    saveas(fig,outfile);
    close(fig);

end


end



function plot_sEEGnal_vs_human(config,bands_info,pow_human_dataset,pow_sEEGnal_dataset)

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
    current_human = pow_human_dataset(:,current_f,:,:);
    current_human = squeeze(nanmean(current_human,2));
    current_sEEGnal = pow_sEEGnal_dataset(:,current_f,:);
    current_sEEGnal = squeeze(nanmean(current_sEEGnal,2));

    % Repmat sEEGnal for each human
    current_sEEGnal = repmat(current_sEEGnal,1,1,size(current_human,3));

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

    % Get the limits for plot purposes
    lims = axis;

    % lsline
    ls = lsline;
    ls.Color = [255 136 136]/255;
    ls.LineWidth = 2;

    % Plot the identity line
    l = line([0 max(lims)],[0 max(lims)],'Color','black','LineStyle','--');

    % Set the limits of the axis to be square
    xlim([0 max(lims)])
    ylim([0 max(lims)])
    axis square

    % Info
    xlabel('Human Expert')
    ylabel('sEEGnal')


end


% Save the figure
outfile = sprintf('%s/pow_sEEGnal_vs_human.svg',config.path.figures);
saveas(fig,outfile);
outfile = sprintf('%s/pow_sEEGnal_vs_human.png',config.path.figures);
saveas(fig,outfile);
close(fig);


end

