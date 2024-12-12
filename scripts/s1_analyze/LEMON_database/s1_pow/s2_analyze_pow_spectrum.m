clear
clc
restoredefaultpath

% Paths
config.path.clean_data = '../../../../databases/LEMON_database/derivatives';
config.path.results = '../../../../results/pow';

if ~exist(config.path.results),mkdir(config.path.results),end

% To define later the pow matrix
complete_channel_labels = {'Fp1', 'Fp2', 'F7', 'F3', 'Fz', 'F4', 'F8', 'FC5', 'FC1', 'FC2', 'FC6', 'T7', 'C3', 'Cz', 'C4', 'T8', 'CP5', 'CP1', 'CP2', 'CP6', 'AFz', 'P7', 'P3', 'Pz', 'P4', 'P8', 'PO9', 'O1', 'Oz', 'O2', 'PO10', 'AF7', 'AF3', 'AF4', 'AF8', 'F5', 'F1', 'F2', 'F6', 'FT7', 'FC3', 'FC4', 'FT8', 'C5', 'C1', 'C2', 'C6', 'TP7', 'CP3', 'CPz', 'CP4', 'TP8', 'P5', 'P1', 'P2', 'P6', 'PO7', 'PO3', 'POz', 'PO4', 'PO8'};

% Get the different testers
testers = {'lemon','sEEGnal'};

% Read the power spectrum
[pow_lemon_dataset_norm,~,channels_lemon_included] = read_pow_dataset(config,'lemon');
[pow_sEEGnal_dataset_norm,f,channels_sEEGnal_included] = read_pow_dataset(config,'sEEGnal');

% Per band and area
% Define the frequency bands
bands_info = struct('name',{'delta', 'theta','alpha','beta','gamma', 'broadband'},...
    'f_limits',{[ 2 4] , [4 8] , [8 12] , [12 20] , [20 45] , [2 45]},...
    'f_limits_index',[],'f_original',f');

% Variable for the stats
results = [];

for iband = 1 : numel(bands_info)
    
    % Obtain the current freqencies of interest
    current_band = bands_info(iband).name;
    current_band_f_limits = bands_info(iband).f_limits;
    current_f = f>current_band_f_limits(1) & f<current_band_f_limits(2);
    current_f = current_f';
    bands_info(iband).f_limits_index = current_f;
    
    % To do the mean and std of the statistics
    stats = [];
    stats.NMSE = nan(numel (complete_channel_labels),size(pow_lemon_dataset_norm,3));
    stats.rho = nan(numel (complete_channel_labels),size(pow_lemon_dataset_norm,3));
    stats.tstat = nan(numel (complete_channel_labels),size(pow_lemon_dataset_norm,3));
    stats.cohen_d = nan(numel (complete_channel_labels),size(pow_lemon_dataset_norm,3));
    stats.p = nan(numel (complete_channel_labels),size(pow_lemon_dataset_norm,3));
    
    % For each sensor
    for ichannel = 1 : numel (complete_channel_labels)
        
        % Get the info for the current comparison
        current_band_pow = nan(sum(current_f),size(pow_lemon_dataset_norm,3),2);
        current_band_pow(:,:,1) = pow_lemon_dataset_norm(ichannel,current_f,:);
        current_band_pow(:,:,2) = pow_sEEGnal_dataset_norm(ichannel,current_f,:);
        
        % For each subject, estimate the NMSE, corr, and t-test comparing
        % the two pre-processing pipelines
        for isubject = 1 : size(current_band_pow,2)
            
            % NMSE
            lemon = current_band_pow(:,isubject,1);
            sEEGnal = current_band_pow(:,isubject,2);
            numerator = sum((lemon - sEEGnal).^2);
            denominator = sum(lemon.^2);
            NMSE = numerator / denominator;
            stats.NMSE(ichannel,isubject) = NMSE;
            
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
    results.stats.(current_band).NMSE = stats.NMSE;
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
outfile = sprintf('%s/pow_results.mat',config.path.results);
save(outfile,'-struct','results');


% Functions
function [pow_dataset_norm,f,channels] = read_pow_dataset(config,dataset_name)

% Load the datset
dataset_path = sprintf('%s/%s/%s_dataset.mat',config.path.clean_data,...
    dataset_name,dataset_name);
dummy = load(dataset_path);

for icurrent = 1 : numel(dummy.dataset)
    
    % Load pow
    current_dataset = dummy.dataset(icurrent);
    pow_file = sprintf('%s/%s',current_dataset.pow.path,...
        current_dataset.pow.file);
    pow = load(pow_file);
    
    % Save f
    f = pow.f;
    
    % Normalize
    current_pow = mean(pow.pow_spectrum,3);
    scaling_factor = 1./nansum(nansum(current_pow,2));
    scaling_factor = repmat(scaling_factor,size(current_pow));
    current_pow_norm = scaling_factor .* current_pow;
    
    % Add to the all matrix
    if icurrent == 1
        pow_dataset_norm = nan(numel(pow.complete_channel_labels),numel(f),numel(dataset));
        channels = struct('channels_included',[],...
            'channels_included_index',[]);
    end
    pow_dataset_norm(1:size(current_pow_norm,1),:,icurrent) = current_pow_norm;
    channels(icurrent).channels_included = pow.channels_included;
    channels(icurrent).channels_included_index = pow.channels_included_index;
    
end

end

