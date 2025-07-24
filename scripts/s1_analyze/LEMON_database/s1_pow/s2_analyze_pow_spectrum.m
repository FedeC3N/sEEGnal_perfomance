%{

Evaluate differnces between sEEGnal and human experts.
 
For each subject, we estimate correlation, NRMSE, t-test sensor by sensor
and save the results.


@author: Fede

%}

clear
clc
restoredefaultpath

% Paths
config.path.clean_data = '../../../../databases/LEMON_database/derivatives';
config.path.results = '../../../../results/LEMON_database/pow';

if ~exist(config.path.results),mkdir(config.path.results),end

% To define later the pow matrix
complete_channel_labels = {'Fp1', 'Fp2', 'F7', 'F3', 'Fz', 'F4', 'F8', 'FC5', 'FC1', 'FC2', 'FC6', 'T7', 'C3', 'Cz', 'C4', 'T8', 'CP5', 'CP1', 'CP2', 'CP6', 'AFz', 'P7', 'P3', 'Pz', 'P4', 'P8', 'PO9', 'O1', 'Oz', 'O2', 'PO10', 'AF7', 'AF3', 'AF4', 'AF8', 'F5', 'F1', 'F2', 'F6', 'FT7', 'FC3', 'FC4', 'FT8', 'C5', 'C1', 'C2', 'C6', 'TP7', 'CP3', 'CPz', 'CP4', 'TP8', 'P5', 'P1', 'P2', 'P6', 'PO7', 'PO3', 'POz', 'PO4', 'PO8'};

% Get the different testers
testers = {'lemon','sEEGnal'};

% Per band and area
% Define the frequency bands
bands_info = struct('name',{'delta', 'theta','alpha','beta','gamma', 'broadband'},...
    'f_limits',{[ 2 4] , [4 8] , [8 12] , [12 20] , [20 45] , [2 45]},...
    'f_limits_index',[],'f_original',[]);

% Variable for the stats
results = [];
stats = []; % dimensions: [channels,recordings,bands]
stats.NRMSE = nan(numel (complete_channel_labels),20,numel(bands_info));
stats.rho = nan(numel (complete_channel_labels),20,numel(bands_info));
stats.tstat = nan(numel (complete_channel_labels),20,numel(bands_info));
stats.cohen_d = nan(numel (complete_channel_labels),20,numel(bands_info));
stats.p = nan(numel (complete_channel_labels),20,numel(bands_info));


% Read the power spectrum
[pow_lemon_dataset_norm,~,channels_human_included] = read_pow_dataset(config,'lemon');
[pow_sEEGnal_dataset_norm,f,channels_sEEGnal_included] = read_pow_dataset(config,'sEEGnal');

for iband = 1 : numel(bands_info)

    % Obtain the current freqencies of interest
    current_band = bands_info(iband).name;
    current_band_f_limits = bands_info(iband).f_limits;
    current_f = f>current_band_f_limits(1) & f<current_band_f_limits(2);
    current_f = current_f';
    bands_info(iband).f_limits_index = current_f;
    bands_info(iband).f_original = f';


    % For each sensor
    for ichannel = 1 : numel (complete_channel_labels)

        % Get the info for the current comparison
        current_band_pow = nan(sum(current_f),size(pow_lemon_dataset_norm,3),2);
        current_band_pow(:,:,1) = pow_lemon_dataset_norm(ichannel,current_f,:);
        current_band_pow(:,:,2) = pow_sEEGnal_dataset_norm(ichannel,current_f,:);

        % For each subject, estimate the NMSE, corr, and t-test comparing
        % the two pre-processing pipelines
        for isubject = 1 : size(current_band_pow,2)

            % NRMSE
            lemon =  current_band_pow(:,isubject,1);
            sEEGnal = current_band_pow(:,isubject,2);
            MSE = nanmean((lemon - sEEGnal).^2);
            RMSE = sqrt(MSE);
            NRMSE = RMSE / 1; % To normalize we use the range, in norm_pov [0 1]
            stats.NRMSE(ichannel,isubject,iband) = NRMSE;

            % corr
            [rho,~] = corrcoef(lemon,sEEGnal,'Rows','complete');
            stats.rho(ichannel,isubject,iband) = rho(1,2);

            % t-test and  Effect Size (Cohen's d)
            diff = lemon - sEEGnal;
            [~,p,~,t_stats] = ttest(diff);
            diff = abs(lemon - sEEGnal);
            d = mean(diff)/std(diff);
            stats.tstat(ichannel,isubject,iband) = t_stats.tstat;
            stats.cohen_d(ichannel,isubject,iband) = d;
            stats.p(ichannel,isubject,iband) = p;

        end

    end


end



% Complete the information
results.stats = stats;
results.dimensions = ['channels x recordings x bands'];
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
    pow_file = sprintf('../../../../%s/%s',current_dataset.pow.path,...
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

