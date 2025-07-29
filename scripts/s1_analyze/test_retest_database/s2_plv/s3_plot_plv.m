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
config.path.clean_data = '../../../../databases/test_retest_database/derivatives';
config.path.results = '../../../../results/test_retest_database/plv';
config.path.figures = '../../../../docs/manuscript/figures/test_retest_database/plv_results';

if ~exist(config.path.figures), mkdir(config.path.figures),end

% Load the results
stats = load(sprintf('%s/plv_results.mat',config.path.results));

% Define the frequency bands
bands_info = struct('name',{'delta', 'theta','alpha','beta','gamma'},...
    'f_limits',{[ 2 4] , [4 8] , [8 12] , [12 20] , [20 45] });
config.bands_info = stats.bands_info;

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
config.bands_info = bands_info;

% Define measures
config.measures = {'NRMSE'};

% Desired tasks
config.tasks = {'EO','EC'};

% Channels
config.complete_channel_labels = stats.complete_channel_labels;

for itask = 1 : numel(config.tasks)

    % Save for plot title purpose
    config.itask = itask;

    % Read the power spectrum
    % Differentiate EC and EO
    plv_dataset = read_plv_dataset(config, config.tasks{itask});
    current_stats = stats.stats;
    current_stats.NRMSE = current_stats.NRMSE(:,:,:,itask);
    current_stats.rho = current_stats.rho(:,:,:,itask);
    current_stats.tstat = current_stats.tstat(:,:,:,itask);
    current_stats.cohen_d = current_stats.cohen_d(:,:,:,itask);
    current_stats.p = current_stats.p(:,:,:,itask);

    %%%%%%%%%%%%%%%
    % PLOTS
    %%%%%%%%%%%%%%%

    % Topoplots of the stats
    plot_stats_in_head(config,current_stats)

    % Plot differences in channels as measures by NRMSE
    plot_violinplot(config,stats.stats,stats.bands_info)

    % Plot differences in channels as measures by corr
    plot_first_vs_second(config,stats.bands_info,plv_dataset)

end

% Aux functions
function plv_dataset = read_plv_dataset(config,desired_task)

% Load the datset
dataset_path = sprintf('%s/sEEGnal/sEEGnal_dataset.mat',...
    config.path.clean_data);
dummy = load(dataset_path);

% Keep only the files of the condition
subjects = {dummy.dataset.sub};
subjects = unique(subjects);

% Create the output
% [channels x channels x bands_info x recordings x subjects]
plv_dataset = nan(numel(config.complete_channel_labels),numel(config.complete_channel_labels),...
    numel(config.bands_info),2,numel(subjects));

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
        plv_file = sprintf('../../../../%s/%s.mat',current_dataset.plv.path,...
            current_dataset.plv.file);
        plv = load(plv_file);

        for iband = 1 : numel(config.bands_info)

            % Add to the all matrix
            plv_dataset(:,:,iband,irecording,isubject) = plv.(config.bands_info(iband).name).plv;

        end


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
    outfile = sprintf('%s/%s_%s_whole_head.svg',...
        config.path.figures,config.tasks{config.itask},...
        current_measure);
    saveas(fig,outfile);
    outfile = sprintf('%s/%s_%s_whole_head.png',...
        config.path.figures,config.tasks{config.itask},...
        current_measure);
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
    outfile = sprintf('%s/%s_plv_%s_violinplot.svg',...
        config.path.figures, config.tasks{config.itask},current_measure);
    saveas(fig,outfile);
    outfile = sprintf('%s/%s_plv_%s_violinplot.png',...
        config.path.figures, config.tasks{config.itask},current_measure);
    saveas(fig,outfile);
    close(fig);

end


end


function plot_first_vs_second(config,bands_info,plv_dataset)

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
    current_first = squeeze(plv_dataset(:,:,iband,1,:));
    % Create the upper matrix mask
    triu_mask = triu(ones(size(current_first(:,:,1))));
    triu_mask = repmat(triu_mask,1,1,size(current_first,3));
    triu_mask = logical(triu_mask);
    current_first = current_first(triu_mask);
    current_second = squeeze(plv_dataset(:,:,iband,2,:));
    current_second = current_second(triu_mask);

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
    xlabel('First recording')
    ylabel('Second recording')

end


% Save the figure
outfile = sprintf('%s/%s_plv_first_vs_second.svg',config.path.figures,...
    config.tasks{config.itask});
saveas(fig,outfile);
close(fig);

end