clear
clc
restoredefaultpath

% Paths
config.path.clean_data = '../../../../databases/AI_Mind_database/derivatives';
config.path.results = '../../../../results/AI_Mind_database/pow';

if ~exist(config.path.results),mkdir(config.path.results),end

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

% Get the different testers (except sEEGnal)
testers = dir(sprintf('%s/*',config.path.clean_data));
testers = testers(3:end-1);

% Define the frequency bands
bands_info = struct('name',{'delta', 'theta','alpha','beta','gamma', 'broadband'},...
    'f_limits',{[ 2 4] , [4 8] , [8 12] , [12 20] , [20 45] , [2 45]},...
    'f_limits_index',[],'f_original',[]);

% Variable for the stats
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
    [pow_human_dataset_norm,~,channels_human_included] = read_pow_dataset(config,current_tester);
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
            current_band_pow = nan(sum(current_f),size(pow_human_dataset_norm,3),2);
            current_band_pow(:,:,1) = pow_human_dataset_norm(ichannel,current_f,:);
            current_band_pow(:,:,2) = pow_sEEGnal_dataset_norm(ichannel,current_f,:);
            
            % For each subject, estimate the NMSE, corr, and t-test comparing
            % the two pre-processing pipelines
            for isubject = 1 : size(current_band_pow,2)
                
                % NRMSE
                human =  current_band_pow(:,isubject,1);
                sEEGnal = current_band_pow(:,isubject,2);
                MSE = nanmean((human - sEEGnal).^2);
                RMSE = sqrt(MSE);
                NRMSE = RMSE / 1; % To normalize we use the range, in norm_pov [0 1]
                stats.NRMSE(ichannel,isubject,iband,itester) = NRMSE;
                
                % corr
                [rho,~] = corrcoef(human,sEEGnal,'Rows','complete');
                stats.rho(ichannel,isubject,iband,itester) = rho(1,2);
                
                % t-test and  Effect Size (Cohen's d)
                diff = human - sEEGnal;
                [~,p,~,t_stats] = ttest(diff);
                diff = abs(human - sEEGnal);
                d = mean(diff)/std(diff);
                stats.tstat(ichannel,isubject,itester) = t_stats.tstat;
                stats.cohen_d(ichannel,isubject,itester) = d;
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
    pow_file = sprintf('../../../../%s/%s.mat',current_dataset.pow.path,...
        current_dataset.pow.file);
    if ~exist(pow_file)
        fake_pow = nan(size(current_pow));
        pow_dataset_norm(1:size(current_pow_norm,1),:,icurrent) = fake_pow;
        channels(icurrent).channels_included = nan;
        channels(icurrent).channels_included_index = nan;
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

