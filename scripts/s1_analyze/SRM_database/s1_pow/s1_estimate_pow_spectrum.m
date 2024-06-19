clear
clc
restoredefaultpath

% Add functions to the path
addpath('../shared/')
addpath('../../../../../../SharedFunctions/fieldtrip-20220626/')
ft_defaults

% Paths
config.path.dataset = '../../../../data/SRM_database/dataset';
config.path.pow = '../../../../data/SRM_database/pow';

if ~exist(config.path.pow), mkdir(config.path.pow), end

% To define later the pow matrix
complete_channel_labels = {'Fp1';'AF7';'AF3';'F1';'F3';'F5';'F7';'FT7';'FC5';'FC3';'FC1';'C1';'C3';'C5';'T7';'TP7';'CP5';'CP3';'CP1';'P1';'P3';'P5';'P7';'P9';'PO7';'PO3';'O1';'Iz';'Oz';'POz';'Pz';'CPz';'Fpz';'Fp2';'AF8';'AF4';'AFz';'Fz';'F2';'F4';'F6';'F8';'FT8';'FC6';'FC4';'FC2';'FCz';'Cz';'C2';'C4';'C6';'T8';'TP8';'CP6';'CP4';'CP2';'P2';'P4';'P6';'P8';'P10';'PO8';'PO4';'O2'};

% Load the whole dataset
load(sprintf('%s/SRM_dataset.mat',config.path.dataset));

for ifile = 1 : numel(dataset)
    
    fprintf('Working on %s\n\n', dataset(ifile).file)
    
    % Process the subject
    [~,data] = evalc( 'process_subject(dataset(ifile))' );
    data_matrix = cat(3,data.trial{:});
    
    % Force to have 4 seconds (4096 samples)
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
        dataset(ifile).database);
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
outfile = sprintf('%s/SRM_dataset.mat',config.path.dataset);
save('-v7.3',outfile,'dataset')

