clear
clc
restoredefaultpath

% Add functions to the path
addpath('../shared/')
addpath('../../../../../../Programacion/SharedFunctions/fieldtrip-20220626/')
ft_defaults

% Paths
config.path.clean_data = '../../../../databases/AI_Mind_database/derivatives';

% To define later the pow matrix
complete_channel_labels = {'Fp1', 'Fp2', 'F7', 'F3', 'Fz', 'F4', 'F8', 'FC5', 'FC1', 'FC2', 'FC6', 'T7', 'C3', 'Cz', 'C4', 'T8', 'CP5', 'CP1', 'CP2', 'CP6', 'AFz', 'P7', 'P3', 'Pz', 'P4', 'P8', 'PO9', 'O1', 'Oz', 'O2', 'PO10', 'AF7', 'AF3', 'AF4', 'AF8', 'F5', 'F1', 'F2', 'F6', 'FT7', 'FC3', 'FC4', 'FT8', 'C5', 'C1', 'C2', 'C6', 'TP7', 'CP3', 'CPz', 'CP4', 'TP8', 'P5', 'P1', 'P2', 'P6', 'PO7', 'PO3', 'POz', 'PO4', 'PO8'};

% Avoid overwrite
config.overwrite = false;

% Get the different testers
testers = dir(sprintf('%s/*',config.path.clean_data));
testers = testers(3:end);

for itester = 1 : numel(testers)
    
    % Load the datset
    current_tester = testers(itester);
    dataset_path = sprintf('%s/%s/%s_dataset.mat',current_tester.folder,...
        current_tester.name,current_tester.name);
    load(dataset_path);
    
    for ifile = 1 : numel(dataset)
        
        % Only .set files
        if ~strcmp(dataset(ifile).file(end-2:end),'set')
            continue
        end
        
        fprintf('Working on %s\n\n', dataset(ifile).file)
        
        % Create the outfile name
        path_pow = sprintf('../../../../databases/LEMON_database/derivatives/%s/pow',...
            current_tester.name);
        if ~exist(path_pow), mkdir(path_pow), end
        outfile_name = sprintf('sub-%s_ses-%s_task-%s_%s_pow',...
            dataset(ifile).sub,dataset(ifile).ses,dataset(ifile).task,...
            dataset(ifile).origin);
        outfile = sprintf('%s/%s.mat',path_pow,outfile_name);
        
        % Check if exist and overwrite
        if ~exist(outfile) || (exist(outfile) && config.overwrite)
            
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
                channels_included_index = false(size(complete_channel_labels));
                for ichannel = 1 : numel(data.label)
                    
                    current_channel_index = find(ismember(complete_channel_labels,...
                        data.label{ichannel}));
                    
                    if ~isempty(current_channel_index)
                        channels_included_index(current_channel_index) = true;
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
            pow.channels_included_index = channels_included_index;
            pow.pow_spectrum = pow_spectrum;
            
            % Save the file
            save(outfile,'-struct','pow')
            
            % Save the metadata in the dataset
            pow = [];
            pow.path = fullfile('databases','LEMON_database','derivatives',...
                current_tester.name, 'pow');;
            pow.file = outfile_name;
            dataset(ifile).pow = pow;
            
        else
            
            fprintf('   Already calculated. Do not overwrite.\n\n')
            
        end
        
        
    end
    
    % Save
    save('-v7.3',dataset_path,'dataset')
    
end