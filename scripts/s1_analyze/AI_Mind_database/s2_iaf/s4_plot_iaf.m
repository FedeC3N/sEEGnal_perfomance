clear
clc
close all
restoredefaultpath

% Add paths
addpath('../shared/')

% Paths
config.path.dataset = '../../../../metadata/AI_Mind_database/dataset';
config.path.stats = '../../../../data/AI_Mind_database/stats';

% Load the whole dataset
load(sprintf('%s/AI_Mind_dataset.mat',config.path.dataset));

% Load the IAF
iaf_eeg_expert_dataset = read_iaf_eeg_expert_dataset(dataset);
iaf_ETL_dataset = read_iaf_dataset(dataset,'etl');

% Measures
measures = {'iaf', 'iaf_amp'};

% Plot difference in channels
plot_diff_channels(iaf_eeg_expert_dataset,iaf_ETL_dataset,measures)

% Aux functions
function iaf_dataset = read_iaf_dataset(dataset, desired_dataset)

% subject of interest
current_dataset_index = ismember({dataset.origin},desired_dataset);
current_dataset = dataset(current_dataset_index);

for icurrent = 1 : numel(current_dataset)
    
    % Load pow
    iaf = load(sprintf('%s/%s',current_dataset(icurrent).iaf.path,...
        current_dataset(icurrent).iaf.file));
    
    % Add to the all matrix
    if icurrent == 1
        iaf_dataset = struct('iaf',[],'iaf_amp',[]);
    end
    iaf_dataset(icurrent).iaf = iaf.iaf;
    iaf_dataset(icurrent).iaf_amp =iaf.iaf_amp;
    
end

end


function iaf_eeg_expert_dataset = read_iaf_eeg_expert_dataset(dataset)

users = {dataset.origin};
users = unique(users(~ismember(users,'etl')));
for iuser = 1 : numel(users)
    
    current_iaf_eeg_expert_dataset = ...
        read_iaf_dataset(dataset,users{iuser});
    
    % Little bit complicated struct
    % For the first eeg_expert, store the result
    if iuser == 1
        iaf_eeg_expert_dataset = struct('iaf',[],'iaf_amp',[]);
        
        for isubject = 1 : numel(current_iaf_eeg_expert_dataset)
            iaf_eeg_expert_dataset(isubject).iaf = current_iaf_eeg_expert_dataset(isubject).iaf;
            iaf_eeg_expert_dataset(isubject).iaf_amp = current_iaf_eeg_expert_dataset(isubject).iaf_amp;
        end
    % For the rest, sum them
    else
        for isubject = 1 : numel(current_iaf_eeg_expert_dataset)
            iaf_eeg_expert_dataset(isubject).iaf = iaf_eeg_expert_dataset(isubject).iaf + current_iaf_eeg_expert_dataset(isubject).iaf;
            iaf_eeg_expert_dataset(isubject).iaf_amp = iaf_eeg_expert_dataset(isubject).iaf_amp + current_iaf_eeg_expert_dataset(isubject).iaf_amp;
        end
    end
    
end

% Finally, obtain tha average dividing by the number of eeg_experts
for isubject = 1 : numel(current_iaf_eeg_expert_dataset)
    iaf_eeg_expert_dataset(isubject).iaf = iaf_eeg_expert_dataset(isubject).iaf/numel(users);
    iaf_eeg_expert_dataset(isubject).iaf_amp = iaf_eeg_expert_dataset(isubject).iaf_amp/numel(users);
end

end


function plot_diff_channels(iaf_eeg_expert_dataset,iaf_ETL_dataset,measures)

figure('WindowState', 'maximized');
hold on

% For each band
for imeasure = 1 : numel(measures)
    
    % Select the current measure
    current_iaf_eeg_expert = [iaf_eeg_expert_dataset(:).(measures{imeasure})];
    current_iaf_ETL = [iaf_ETL_dataset(:).(measures{imeasure})];
    
    % Get the difference
    iaf_diff_vector = current_iaf_eeg_expert - current_iaf_ETL;
    x_vector = imeasure * ones(numel(iaf_diff_vector),1);
    
    
    % Plot
    sw = swarmchart(x_vector,iaf_diff_vector,'filled','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5);
    bx = boxchart(x_vector,iaf_diff_vector,'BoxFaceColor',sw.CData ./ 1.2,'WhiskerLineColor',sw.CData ./ 1.2,...
        'MarkerStyle','none','BoxWidth',sw.XJitterWidth);
    
    
end

xlim([0 numel(measures) + 1])
xticks(1:numel(measures))
xticklabels(measures)
title('Average IAF difference')
set(gca,'TickLabelInterpreter','none')



end


function plot_corr_channels(plv_eeg_expert_dataset,plv_ETL_dataset,bands_info, channels)

figure('WindowState', 'maximized');
hold on

% For each band
for iband = 1 : numel(bands_info)
    
    % Select the band data
    current_plv_global_eeg_expert = squeeze(plv_eeg_expert_dataset(:,:,iband,:));
    current_plv_global_ETL = squeeze(plv_ETL_dataset(:,:,iband,:));
    plv_vector_eeg_expert = reshape(current_plv_global_eeg_expert,[],size(current_plv_global_eeg_expert,3));
    plv_vector_ETL = reshape(current_plv_global_ETL,[],size(current_plv_global_ETL,3));
    
    % Get the upper matrix
    dummy = ones(size(current_plv_global_ETL(:,:,1)));
    upper_index = abs(triu(dummy,1)) > 0 ;
    upper_index = upper_index(:);
    plv_vector_eeg_expert =  plv_vector_eeg_expert(upper_index,:);
    plv_vector_ETL = plv_vector_ETL(upper_index,:);
    
    % Estimate correlation for each channel
    plv_vector_eeg_expert = plv_vector_eeg_expert';
    plv_vector_ETL = plv_vector_ETL';
    rho = nan(1,size(plv_vector_ETL,2));
    x_vector = iband * ones(1,size(plv_vector_ETL,2));
    for ichannel = 1 : size(plv_vector_ETL,2)
        
        % Remove nans
        current_channel_eeg_expert = plv_vector_eeg_expert(:,ichannel);
        current_channel_ETL = plv_vector_ETL(:,ichannel);
        nan_index = isnan(current_channel_eeg_expert) | isnan(current_channel_ETL);
        current_channel_eeg_expert = current_channel_eeg_expert(~nan_index);
        current_channel_ETL = current_channel_ETL(~nan_index);
        
        
        rho(ichannel) = corr(current_channel_eeg_expert,current_channel_ETL);
        
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


function plot_corr_plv_subjects(plv_eeg_expert_dataset,plv_ETL_dataset,bands_info, channels_eeg_expert,channels_ETL)

figure('WindowState', 'maximized');
hold on

% For each band
for iband = 1 : numel(bands_info)
    
    % Create the rho empty and the plot axis
    rho = nan(1,size(plv_eeg_expert_dataset,4));
    x_vector = iband * ones(1,size(rho,2));
    for irho = 1 :numel(rho)
        
        % Get the current plv matrix
        current_plv_eeg_expert = plv_eeg_expert_dataset(:,:,iband,irho);
        current_plv_ETL = plv_ETL_dataset(:,:,iband,irho);
        
        % Create an index of valid channels
        valid = channels_eeg_expert(irho).channels_included_index & channels_ETL(irho).channels_included_index;
        
        % Remove badchannels
        current_plv_eeg_expert = current_plv_eeg_expert(valid,valid);
        current_plv_ETL = current_plv_ETL(valid,valid);
        
        % Get the upper matrix in vector format
        dummy = ones(size(current_plv_ETL));
        upper_index = abs(triu(dummy,1)) > 0 ;
        current_plv_eeg_expert = current_plv_eeg_expert(upper_index);
        current_plv_ETL = current_plv_ETL(upper_index);
        
        % Correlation
        rho(irho) = corr(current_plv_eeg_expert,current_plv_ETL);
        
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


