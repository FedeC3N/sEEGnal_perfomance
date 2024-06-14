% PLVs needs the channels information and transform the matrix into 
% 64 x 64
clear
clc
restoredefaultpath

% Add functions to the path
addpath('../../../../SharedFunctions/fieldtrip-20220626/')
ft_defaults

% Paths
config.path.dataset = '../../data/SRM_database/dataset';
config.path.plv = '../../data/SRM_database/plv_bad';

% To define later the PLV matrix
complete_channel_labels = {'Fp1';'AF7';'AF3';'F1';'F3';'F5';'F7';'FT7';'FC5';'FC3';'FC1';'C1';'C3';'C5';'T7';'TP7';'CP5';'CP3';'CP1';'P1';'P3';'P5';'P7';'P9';'PO7';'PO3';'O1';'Iz';'Oz';'POz';'Pz';'CPz';'Fpz';'Fp2';'AF8';'AF4';'AFz';'Fz';'F2';'F4';'F6';'F8';'FT8';'FC6';'FC4';'FC2';'FCz';'Cz';'C2';'C4';'C6';'T8';'TP8';'CP6';'CP4';'CP2';'P2';'P4';'P6';'P8';'P10';'PO8';'PO4';'O2'};

% Load the whole dataset
load(sprintf('%s/SRM_dataset.mat',config.path.dataset));

for idataset = 1 : numel(dataset)
    
    % Load the current file
    current_dataset = dataset(idataset);
    current_plv_file = sprintf('%s/%s',config.path.plv,current_dataset.plv.file);
    plv = load(current_plv_file);
    
    fprintf('Working on %s\n\n', current_dataset.file)
    
    % Remove broadband
    plv.bands_info(end) = [];
    plv = rmfield(plv,'broadband');
    
    % Save
    save(current_plv_file,'-struct','plv')
    
    continue
    
    % Load the subject data
    current_dataset = dataset(idataset);
    current_dataset.path = current_dataset.path(7:end);
    data = process_subject(current_dataset);
    
    % Load the current file
    current_plv_file = sprintf('%s/%s',config.path.plv,current_dataset.plv.file);
    plv = load(current_plv_file);
    
    % Create an index of the channels position
    channels_included = data.label;
    channels_included_index = ismember(complete_channel_labels,channels_included);
    
    % Save the information
    plv.channels_included = channels_included;
    plv.channels_included_index = channels_included_index;
    
    % Correct each band
    for iband = 1 : numel(plv.bands_info)
        
        new_plv = nan(64,64);
        old_plv = plv.(plv.bands_info(iband).name).plv;
        
        % Avoid the non-included channels
        new_plv(channels_included_index,channels_included_index) = old_plv;
        
        % Save it
        plv.(plv.bands_info(iband).name).plv = new_plv;
        
        
    end
    
    % Save it
    save(current_plv_file,'-struct','plv');
    
    
    
    
end


% Aux functions
function data = process_subject(current_dataset)
    
    dummy_complete_file = sprintf('%s/%s',current_dataset.path,...
        current_dataset.file);
    
    % Load the file
    cfg = [];
    cfg.dataset = dummy_complete_file;
    data = ft_preprocessing(cfg);
    
    % Process the data
    switch current_dataset.database
        
        case 'SRM_database'
            
            % Read the badchannels
            badchannel_file = dir(sprintf('%s/*channels.tsv',...
                current_dataset.path));
            badchannel_file = sprintf('%s/%s',current_dataset.path,...
                badchannel_file.name);
            
            % Fields of interest: name and status
            tsv_content = tdfread(badchannel_file);
            
            % Find the badchannels
            bad_index = ismember(cellstr(tsv_content.status),'bad');
            badchannels = cellstr(tsv_content.name);
            badchannels = badchannels(bad_index);
            
            % Remove the badchannels
            cfg = [];
            badchannels          = strcat('-',badchannels);
            desired_channel      = cat(2,{'eeg'},badchannels');
            cfg.channel          = desired_channel;
            cfg.precision        = 'single';
            cfg.feedback         = 'no';

            data            = ft_preprocessing ( cfg, data );
            
            
    end


end
