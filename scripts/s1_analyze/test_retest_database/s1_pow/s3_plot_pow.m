%{

Different plots to investigate the results.


@author: Fede

%}

clear
clc
close all
restoredefaultpath

% Paths
config.path.clean_data = '../../../../databases/test_retest_database/derivatives';
config.path.results = '../../../../results/test_retest_database/pow';
config.path.figures = '../../../../docs/manuscript/figures/test_retest_database/pow_results';

% Add paths
addpath('../shared/')

if ~exist(config.path.figures), mkdir(config.path.figures), end

% Load the results
stats = load(sprintf('%s/pow_results.mat',config.path.results));

% Get the different testers (except sEEGnal)
testers = dir(sprintf('%s/*',config.path.clean_data));
testers = testers(3);
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
areas_info(7).channel = stats.complete_channel_labels;

% Define the frequency bands
config.bands = {'delta', 'theta','alpha','beta','gamma', 'broadband'};
config.bands_info = stats.bands_info;

% Define measures
% config.measures = {'NRMSE', 'rho', 'tstat'};
config.measures = {'NRMSE'};

% Desired tasks
config.tasks = {'EO'};

for itask = 1 : numel(config.tasks)

    % Save for plot title purpose
    config.itask = itask;

    % Read the power spectrum
    % Differentiate EC and EO
    [pow_dataset_norm,f,~] = read_pow_dataset(config, config.tasks{itask});
    current_stats = stats.stats;
    current_stats.NRMSE = current_stats.NRMSE(:,:,:,itask);
    current_stats.rho = current_stats.rho(:,:,:,itask);
    current_stats.tstat = current_stats.tstat(:,:,:,itask);
    current_stats.cohen_d = current_stats.cohen_d(:,:,:,itask);
    current_stats.p = current_stats.p(:,:,:,itask);

    %%%%%%%%%%%%%%%
    % PLOTS
    %%%%%%%%%%%%%%%

    % corr and NRMSE in head
    plot_stats_in_head(config,current_stats)

    % Normalized pow spectrum for all areas in broadband
    plot_pow_spectrum_norm_one_subject(config,pow_dataset_norm);

    % Violinplots
    plot_violinplot(config,stats.stats,stats.bands_info)

    % corr
    plot_first_vs_second(config,stats.bands_info,pow_dataset_norm)

end

% Functions
function [pow_dataset_norm,f,channels] = read_pow_dataset(config,desired_task)

% Load the datset
dataset_path = sprintf('%s/sEEGnal/sEEGnal_dataset.mat',...
    config.path.clean_data);
dummy = load(dataset_path);

% Keep only the files of the condition
subjects = {dummy.dataset.sub};
subjects = unique(subjects);

% Create the output
% [channels x freqs x recordings x subjects]
pow_dataset_norm = nan(numel(config.complete_channel_labels),...
    172,2,numel(subjects));
% [channels x recordings x subjects]
channels = false(numel(config.complete_channel_labels),2,numel(subjects));


for isubject = 1 : numel(subjects)

    % Find the subject and 2 recordings of the current task
    dummy_subjects = {dummy.dataset.sub};
    dummy_tasks    = {dummy.dataset.task};
    dummy_tasks    = cellfun(@(x) x(2:end),dummy_tasks,'UniformOutput',false);
    subject_mask = ismember(dummy_subjects,subjects{isubject});
    task_mask = ismember(dummy_tasks,desired_task);
    dummy_mask = subject_mask & task_mask;
    dummy_index = find(dummy_mask);

    if numel(dummy_index) ~= 2
        continue
    end

    for irecording = 1 : numel(dummy_index)

        % Load pow
        current_dataset = dummy.dataset(dummy_index(irecording));
        pow_file = sprintf('../../../../%s/%s.mat',current_dataset.pow.path,...
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
        pow_dataset_norm(:,:,irecording,isubject) = current_pow_norm;
        channels(:,irecording,isubject) = pow.channels_included_mask;


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
        current_stats = stats.(current_measure)(:,:,iband);
        color_elec = nanmean(current_stats,2);

        % Scatter the values of interest
        scatter(pos_elec(:,1),pos_elec(:,2),size_elec,color_elec,'filled');
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
    sgtitle(sprintf('Task %s - Measure %s',config.tasks{config.itask},...
        current_measure))

    % Save the figure
    outfile = sprintf('%s/%s_%s_%s_whole_head_head_measures.svg',...
        config.path.figures,config.tasks{config.itask},...
        current_measure, current_band);
    saveas(fig,outfile);
    outfile = sprintf('%s/%s_%s_%s_whole_head_head_measures.png',...
        config.path.figures,config.tasks{config.itask},...
        current_measure,current_band);
    saveas(fig,outfile);
    close(fig);

end

end



function plot_pow_spectrum_norm_one_subject(config,pow_dataset_norm)


for isubject = 1 : size(pow_dataset_norm,4)

    % Estimate mean error
    fig = figure('WindowState','maximized');
    hold on
    f_of_interest = config.bands_info(6).f_original;

    % Estimate the powers
    current_pow_first = squeeze(pow_dataset_norm(:,:,1,isubject));
    first_average = nanmean(current_pow_first,1);
    current_pow_second = squeeze(pow_dataset_norm(:,:,2,isubject));
    second_average = nanmean(current_pow_second,1);

    % Plot channels first
    plot(f_of_interest,current_pow_first,...
        'Color',[0.6980    0.8902    0.8902],'LineWidth',0.5)
    plot(f_of_interest,current_pow_second,...
        'Color',[1.0000    0.7725    0.6902],'LineWidth',0.5)

    % Plot mean
    plot(f_of_interest,first_average,...
        'Color',[0 0.5 0.5],'LineWidth',3)
    plot(f_of_interest,second_average,...
        'Color',[1 0.498 0.314],'LineWidth',3)

    % Enhance the plot
    lines = findobj(gca, 'Type', 'Line');
    legend(lines(1:2),{'First recording','Second recording'});
    title(sprintf('Task %s - Power in BroadBand in whole head',config.tasks{config.itask}))

    % Save the figure
    outfile = sprintf('%s/%s_%i_pow_spectrum.svg',...
        config.path.figures, config.tasks{config.itask},isubject);
    saveas(fig,outfile);
    outfile = sprintf('%s/%s_%i_pow_spectrum.png',...
        config.path.figures, config.tasks{config.itask}, isubject);
    saveas(fig,outfile);
    close(fig);

end

end



function plot_violinplot(config,stats,bands_info)

% Remove the broadband
bands_info = bands_info(1:end-1);

for imeasure = 1 : numel(config.measures)

    fig = figure('WindowState', 'maximized');
    hold on

    current_measure = config.measures{imeasure};
    current_stats = stats.(current_measure);

    % For each band
    for iband = 1 : numel(bands_info)

        current_band = bands_info(iband).name;
        current_stats_band = current_stats(:,:,iband,config.itask);
        current_stats_band = current_stats_band(:);

        % X axis for plot
        x_vector = iband * ones(numel(current_stats_band),1);

        % Plot
        sw = swarmchart(x_vector,current_stats_band,'filled','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5);
        bx = boxchart(x_vector,current_stats_band,'BoxFaceColor',sw.CData ./ 1.2,'WhiskerLineColor',sw.CData ./ 1.2,...
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
    outfile = sprintf('%s/%s_pow_%s_violinplot.svg',...
        config.path.figures, config.tasks{config.itask},current_measure);
    saveas(fig,outfile);
    outfile = sprintf('%s/%s_pow_%s_violinplot.png',...
        config.path.figures, config.tasks{config.itask},current_measure);
    saveas(fig,outfile);
    close(fig);

end


end



function plot_first_vs_second(config,bands_info,pow_dataset_norm)

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

    current_first = squeeze(pow_dataset_norm(:,current_f,1,:));
    current_first = squeeze(nanmean(current_first,2));
    current_first = current_first(:);
    current_second = squeeze(pow_dataset_norm(:,current_f,2,:));
    current_second = squeeze(nanmean(current_second,2));
    current_second = current_second(:);

    % Remove the nans
    nans_mask = isnan(current_first) | isnan(current_second);
    current_first = current_first(~nans_mask);
    current_second = current_second(~nans_mask);

    % Plot
    subplot(2,3,iband)
    s = scatter(current_first,current_second,'filled');
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
    xlabel('First Recording')
    ylabel('Second recording')


end
sgtitle(config.tasks{config.itask})


% Save the figure
outfile = sprintf('%s/%s_pow_first_vs_second.svg',...
    config.path.figures,config.tasks{config.itask});
saveas(fig,outfile);
outfile = sprintf('%s/%s_pow_first_vs_second.png',...
    config.path.figures,config.tasks{config.itask});
saveas(fig,outfile);
close(fig);


end

