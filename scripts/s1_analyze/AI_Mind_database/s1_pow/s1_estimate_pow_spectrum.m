clear
clc
restoredefaultpath

% Add functions to the path
addpath('../shared/')
addpath('../../../../../../SharedFunctions/fieldtrip-20220626/')
ft_defaults

% Paths
config.path.dataset = '../../../../data/AI_Mind_database/dataset';
config.path.pow = '../../../../data/AI_Mind_database/pow';

if ~exist(config.path.pow), mkdir(config.path.pow), end

% To define later the pow matrix
complete_channel_labels = {'Fp1', 'Fpz', 'Fp2', 'F7', 'F3', 'Fz', 'F4', 'F8', 'FC5', 'FC1', 'FC2', 'FC6', 'M1', 'T7', 'C3', 'Cz', 'C4', 'T8', 'M2', 'CP5', 'CP1', 'CP2', 'CP6', 'P7', 'P3', 'Pz', 'P4', 'P8', 'POz', 'O1', 'O2', 'AF7', 'AF3', 'AF4', 'AF8', 'F5', 'F1', 'F2', 'F6', 'FC3', 'FCz', 'FC4', 'C5', 'C1', 'C2', 'C6', 'CP3', 'CP4', 'P5', 'P1', 'P2', 'P6', 'F9', 'PO3', 'PO4', 'F10', 'FT7', 'FT8', 'TP7', 'TP8', 'PO7', 'PO8', 'FT9', 'FT10', 'TPP9h', 'TPP10h', 'PO9', 'PO10', 'P9', 'P10', 'AFF1', 'AFz', 'AFF2', 'FFC5h', 'FFC3h', 'FFC4h', 'FFC6h', 'FCC5h', 'FCC3h', 'FCC4h', 'FCC6h', 'CCP5h', 'CCP3h', 'CCP4h', 'CCP6h', 'CPP5h', 'CPP3h', 'CPP4h', 'CPP6h', 'PPO1', 'PPO2', 'I1', 'Iz', 'I2', 'AFp3h', 'AFp4h', 'AFF5h', 'AFF6h', 'FFT7h', 'FFC1h', 'FFC2h', 'FFT8h', 'FTT9h', 'FTT7h', 'FCC1h', 'FCC2h', 'FTT8h', 'FTT10h', 'TTP7h', 'CCP1h', 'CCP2h', 'TTP8h', 'TPP7h', 'CPP1h', 'CPP2h', 'TPP8h', 'PPO9h', 'PPO5h', 'PPO6h', 'PPO10h', 'POO9h', 'POO3h', 'POO4h', 'POO10h', 'OI1h', 'OI2h'};

% Load the whole dataset
load(sprintf('%s/AI_Mind_dataset.mat',config.path.dataset));

for ifile = 1 : numel(dataset)
    
    fprintf('Working on %s\n\n', dataset(ifile).file)
    
    % Process the subject
    [~,data] = evalc( 'process_subject(dataset(ifile))' );
    data_matrix = cat(3,data.trial{:});
    
    % Force to have 4 seconds
    n = 1:4*data.hdr.Fs;
    data_matrix = data_matrix(:,n,:);

    % Define the parameters for power spectrum estimation
    windows_size = size(data_matrix,2);
    noverlap = windows_size/2;
    nfft = size(data_matrix,2)-1;
    fs = data.hdr.Fs;
    pow_spectrum = nan(numel(complete_channel_labels),windows_size/2,size(data_matrix,3));
    
    % Estimate the pow spectrum in each trial
    for itrial = 1 : size(data_matrix,3)
       
        % Matlab works column-wise
        current_trial = data_matrix(:,:,itrial)';
        
        % Estimate the power spectrum (one-sided)
        [current_pow_spectrum,f] = pwelch(current_trial,windows_size,noverlap,nfft,fs,'onesided','power');
        
        % Save in the correct position according to complete channels
        for ichannel = 1 : numel(data.label)
            
            current_channel_index = find(ismember(complete_channel_labels,...
                data.label{ichannel}));
            
            if ~isempty(current_channel_index)
                pow_spectrum(current_channel_index,:,itrial) = current_pow_spectrum(:,ichannel)';
            end
        end
    end
    
    % Keep only frequencies between 2-45 Hz
    f_index = f>2 & f<45;
    pow_spectrum = pow_spectrum(:,f_index,:);
    f = f(f_index);
    
    % Add to the dataset struct
    pow = [];    
    pow.f = f;
    pow.complete_channel_labels = complete_channel_labels;
    pow.channels_included = data.label;
    pow.pow_spectrum = pow_spectrum;
    
    % Save the file
    outfile = sprintf('%s/sub-%s_ses-%s_task-%s_%s_pow', config.path.pow,...
        dataset(ifile).sub,dataset(ifile).ses,dataset(ifile).task,...
        dataset(ifile).origin);
    save(outfile,'-struct','pow')
    
    % Save the metadata in the dataset
    pow = [];
    pow.path = config.path.pow;
    pow.file = sprintf('sub-%s_ses-%s_task-%s_%s_pow',...
        dataset(ifile).sub,dataset(ifile).ses,dataset(ifile).task,...
        dataset(ifile).origin);
    dataset(ifile).pow = pow;
    
end

% Save
outfile = sprintf('%s/AI_Mind_dataset.mat',config.path.dataset);
save('-v7.3',outfile,'dataset')

