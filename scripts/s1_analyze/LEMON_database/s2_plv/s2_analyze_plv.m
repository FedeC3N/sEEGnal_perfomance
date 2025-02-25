clear
clc
restoredefaultpath

% Add paths
addpath('../shared/')

% Paths
config.path.clean_data = '../../../../databases/LEMON_database/derivatives';
config.path.results = '../../../../results/LEMON_database/plv';
if ~exist(config.path.results), mkdir(config.path.results),end

% Get the different testers
testers = {'lemon','sEEGnal'};

% Read the power spectrum
[plv_lemon_dataset,~,channels_lemon_included] = read_plv_dataset(config,'lemon');
[plv_sEEGnal_dataset,f,channels_sEEGnal_included] = read_plv_dataset(config,'sEEGnal');

% Areas
complete_channel_labels = {'Fp1', 'Fp2', 'F7', 'F3', 'Fz', 'F4', 'F8', 'FC5', 'FC1', 'FC2', 'FC6', 'T7', 'C3', 'Cz', 'C4', 'T8', 'CP5', 'CP1', 'CP2', 'CP6', 'AFz', 'P7', 'P3', 'Pz', 'P4', 'P8', 'PO9', 'O1', 'Oz', 'O2', 'PO10', 'AF7', 'AF3', 'AF4', 'AF8', 'F5', 'F1', 'F2', 'F6', 'FT7', 'FC3', 'FC4', 'FT8', 'C5', 'C1', 'C2', 'C6', 'TP7', 'CP3', 'CPz', 'CP4', 'TP8', 'P5', 'P1', 'P2', 'P6', 'PO7', 'PO3', 'POz', 'PO4', 'PO8'};

% Define the frequency bands
bands_info = struct('name',{'delta', 'theta','alpha','beta','gamma'},...
    'f_limits',{[ 2 4] , [4 8] , [8 12] , [12 20] , [20 45] });

% Stats
% Output struct
results = [];

for iband = 1 : numel(bands_info)
    
    % For saving purposes
    current_band = bands_info(iband).name;
    
    % To do the mean and std of the statistics
    stats = [];
    stats.NMSE = nan(numel (complete_channel_labels),size(plv_lemon_dataset,4));
    stats.rho = nan(numel (complete_channel_labels),size(plv_lemon_dataset,4));
    stats.tstat = nan(numel (complete_channel_labels),size(plv_lemon_dataset,4));
    stats.cohen_d = nan(numel (complete_channel_labels),size(plv_lemon_dataset,4));
    stats.p = nan(numel (complete_channel_labels),size(plv_lemon_dataset,4));
    
    for ichannel = 1 : numel(complete_channel_labels)
        
        % Get the current band matrix
        current_band_plv_lemon = squeeze(plv_lemon_dataset(:,:,iband,:));
        current_band_plv_sEEGnal = squeeze(plv_sEEGnal_dataset(:,:,iband,:));
        
        % Get the current PLV for a channel with the rest of the channels
        % Just interested in interconnectivity, not intra-connectivity
        current_band_plv_lemon = squeeze(current_band_plv_lemon(:,ichannel,:));
        current_band_plv_lemon(ichannel,:) = [];
        current_band_plv_sEEGnal = squeeze(current_band_plv_sEEGnal(:,ichannel,:));
        current_band_plv_sEEGnal(ichannel,:) = [];
        
        % For each subject, estimate the NMSE, corr, and t-test comparing
        % the two pre-processing pipelines
        for isubject = 1 : size(current_band_plv_lemon,2)
            
            % NMSE
            lemon = current_band_plv_lemon(:,isubject);
            sEEGnal = current_band_plv_sEEGnal(:,isubject);
            MSE = nanmean((lemon - sEEGnal).^2);
            RMSE = sqrt(MSE);
            NRMSE = RMSE / 2; % To normalize we use the range, in PLV [-1 1]
            stats.NRMSE(ichannel,isubject) = NRMSE;
            
            % corr
            [rho,~] = corrcoef(lemon,sEEGnal,'Rows','complete');
            stats.rho(ichannel,isubject) = rho(1,2);
            
            % t-test and  Effect Size (Cohen's d)
            diff = lemon - sEEGnal;
            [~,p,~,t_stats] = ttest(diff);
            diff = abs(lemon - sEEGnal);
            d = mean(diff)/std(diff);
            stats.tstat(ichannel,isubject) = t_stats.tstat;
            stats.cohen_d(ichannel,isubject) = d;
            stats.p(ichannel,isubject) = p;
            
        end
        
        
    end
    
    % Save
    results.stats.(current_band).NRMSE = stats.NRMSE;
    results.stats.(current_band).rho = stats.rho;
    results.stats.(current_band).tstat = stats.tstat;
    results.stats.(current_band).cohen_d = stats.cohen_d;
    results.stats.(current_band).p = stats.p;
    
end

% Complete the information
results.testers = testers;
results.complete_channel_labels = complete_channel_labels;
results.bands_info = bands_info;

% Save the file
outfile = sprintf('%s/plv_results.mat',config.path.results);
save(outfile,'-struct','results');


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

