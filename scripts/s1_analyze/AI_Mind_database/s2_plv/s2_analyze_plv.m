clear
clc
restoredefaultpath

% Add paths
addpath('../shared/')

% Paths
config.path.clean_data = '../../../../databases/AI_Mind_database/derivatives';
config.path.results = '../../../../results/AI_Mind_database/plv';
if ~exist(config.path.results), mkdir(config.path.results),end

% Get the different testers (except sEEGnal)
testers = dir(sprintf('%s/*',config.path.clean_data));
testers = testers(3:end-1);

% To define later the pow matrix
complete_channel_labels = {'Fp1', 'Fpz', 'Fp2', 'F7', 'F3', 'Fz', 'F4', 'F8', 'FC5', 'FC1', 'FC2', 'FC6',...
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

% Define the frequency bands
bands_info = struct('name',{'delta', 'theta','alpha','beta','gamma'},...
    'f_limits',{[ 2 4] , [4 8] , [8 12] , [12 20] , [20 45] });

% Stats
% Output struct
results = [];
stats = []; % dimensions: [channels,recordings,bands,testers]
stats.NRMSE = nan(numel (complete_channel_labels),20,numel(bands_info),numel(testers));
stats.rho = nan(numel (complete_channel_labels),20,numel(bands_info),numel(testers));
stats.tstat = nan(numel (complete_channel_labels),20,numel(bands_info),numel(testers));
stats.cohen_d = nan(numel (complete_channel_labels),20,numel(bands_info),numel(testers));
stats.p = nan(numel (complete_channel_labels),20,numel(bands_info),numel(testers));

for itester = 1 : numel(testers)
    
    % Read the power spectrum
    current_tester = testers(itester).name;
    [plv_human_dataset,~,channels_lemon_included] = read_plv_dataset(config,current_tester);
    [plv_sEEGnal_dataset,f,channels_sEEGnal_included] = read_plv_dataset(config,'sEEGnal');
    
    for iband = 1 : numel(bands_info)
        
        % For saving purposes
        current_band = bands_info(iband).name;
        
        for ichannel = 1 : numel(complete_channel_labels)
            
            % Get the current band matrix
            current_band_plv_human = squeeze(plv_human_dataset(:,:,iband,:));
            current_band_plv_sEEGnal = squeeze(plv_sEEGnal_dataset(:,:,iband,:));
            
            % Get the current PLV for a channel with the rest of the channels
            % Just interested in interconnectivity, not intra-connectivity
            current_band_plv_human = squeeze(current_band_plv_human(:,ichannel,:));
            current_band_plv_human(ichannel,:) = [];
            current_band_plv_sEEGnal = squeeze(current_band_plv_sEEGnal(:,ichannel,:));
            current_band_plv_sEEGnal(ichannel,:) = [];
            
            % For each subject, estimate the NMSE, corr, and t-test comparing
            % the two pre-processing pipelines
            for isubject = 1 : size(current_band_plv_human,2)
                
                % NMSE
                human = current_band_plv_human(:,isubject);
                sEEGnal = current_band_plv_sEEGnal(:,isubject);
                MSE = nanmean((human - sEEGnal).^2);
                RMSE = sqrt(MSE);
                NRMSE = RMSE / 2; % To normalize we use the range, in PLV [-1 1]
                stats.NRMSE(ichannel,isubject,iband,itester) = NRMSE;
                
                % corr
                [rho,~] = corrcoef(human,sEEGnal,'Rows','complete');
                stats.rho(ichannel,isubject,iband,itester) = rho(1,2);
                
                % t-test and  Effect Size (Cohen's d)
                diff = human - sEEGnal;
                [~,p,~,t_stats] = ttest(diff);
                diff = abs(human - sEEGnal);
                d = mean(diff)/std(diff);
                stats.tstat(ichannel,isubject,iband,itester) = t_stats.tstat;
                stats.cohen_d(ichannel,isubject,iband,itester) = d;
                stats.p(ichannel,isubject,iband,itester) = p;
                
            end
            
            
        end
       
    end
    
end

% Complete the information
results.stats = stats;
results.dimensions = ['channels','recordings','bands','testers'];
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

