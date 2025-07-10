%{

Evaluate consistency intra-session.
 
For each subject, we estimate correlation, NRMSE, t-test sensor by sensor
and save the results.


@author: Fede

%}

clear
clc
restoredefaultpath

% Paths
config.path.clean_data = '../../../../databases/test_retest_database/derivatives';
config.path.results = '../../../../results/test_retest_database/pow';

if ~exist(config.path.results),mkdir(config.path.results),end

% To define later the pow matrix
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

% Get the different testers (except sEEGnal)
testers = dir(sprintf('%s/*',config.path.clean_data));
testers = testers(3);

% Define the frequency bands
bands_info = struct('name',{'delta', 'theta','alpha','beta','gamma', 'broadband'},...
    'f_limits',{[ 2 4] , [4 8] , [8 12] , [12 20] , [20 45] , [2 45]},...
    'f_limits_index',[],'f_original',[]);

% Desired tasks
desired_tasks = {'EO', 'EC'};

% Variable for the stats
results = [];
stats = []; % dimensions: [channels,subjects,bands,task]
stats.NRMSE = nan(numel (config.complete_channel_labels),10,numel(bands_info),numel(desired_tasks));
stats.rho = nan(numel (config.complete_channel_labels),10,numel(bands_info),numel(desired_tasks));
stats.tstat = nan(numel (config.complete_channel_labels),10,numel(bands_info),numel(desired_tasks));
stats.cohen_d = nan(numel (config.complete_channel_labels),10,numel(bands_info),numel(desired_tasks));
stats.p = nan(numel (config.complete_channel_labels),10,numel(bands_info),numel(desired_tasks));

for itask = 1 : numel(desired_tasks)

    current_task = desired_tasks{itask};

    % Differentiate EC and EO
    [pow_dataset_norm,f,~] = read_pow_dataset(config, current_task);

    for iband = 1 : numel(bands_info)

        % Obtain the current freqencies of interest
        current_band = bands_info(iband).name;
        current_band_f_limits = bands_info(iband).f_limits;
        current_f = f>current_band_f_limits(1) & f<current_band_f_limits(2);
        current_f = current_f';
        bands_info(iband).f_limits_index = current_f;
        bands_info(iband).f_original = f';

        % For each sensor
        for ichannel = 1 : numel (config.complete_channel_labels)

            % Get the info for the current comparison
            current_band_pow = nan(sum(current_f),size(pow_dataset_norm,4),2);
            current_band_pow(:,:,1) = squeeze(pow_dataset_norm(ichannel,current_f,1,:));
            current_band_pow(:,:,2) = squeeze(pow_dataset_norm(ichannel,current_f,2,:));

            % For each subject, estimate the NMSE, corr, and t-test comparing
            % the two pre-processing pipelines
            for isubject = 1 : size(current_band_pow,2)

                % NRMSE
                first =  current_band_pow(:,isubject,1);
                second = current_band_pow(:,isubject,2);
                MSE = nanmean((first - second).^2);
                RMSE = sqrt(MSE);
                NRMSE = RMSE / 1; % To normalize we use the range, in norm_pov [0 1]
                stats.NRMSE(ichannel,isubject,iband,itask) = NRMSE;

                % corr
                [rho,~] = corrcoef(first,second,'Rows','complete');
                stats.rho(ichannel,isubject,iband,itask) = rho(1,2);

                % t-test and  Effect Size (Cohen's d)
                diff = first - second;
                [~,p,~,t_stats] = ttest(diff);
                diff = abs(first - second);
                d = mean(diff)/std(diff);
                stats.tstat(ichannel,isubject,itask) = t_stats.tstat;
                stats.cohen_d(ichannel,isubject,itask) = d;
                stats.p(ichannel,isubject,iband,itask) = p;

            end

        end


    end

end

% Complete the information
results.stats = stats;
results.dimensions = ['channels','recordings','bands','task'];
results.testers = testers;
results.config.complete_channel_labels = config.complete_channel_labels;
results.bands_info = bands_info;

% Save the file
outfile = sprintf('%s/pow_results.mat',config.path.results);
save(outfile,'-struct','results');


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

