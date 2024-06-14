% Get the metainformation of each expert and export the clean EEG
% recordings to "cleaned/"
clear
clc
restoredefaultpath

% Paths
config.path.users = '../../data/AI_Mind_database/EEG_experts_input';
config.path.curated = '../../data/AI_Mind_database/curated';
config.path.cleaned = '../../data/AI_Mind_database/cleaned';

% Add fieldtrip and functions
addpath ( '../../../../SharedFunctions/functions/');
addpath('../../../../SharedFunctions/fieldtrip-20220626');
ft_defaults

% Read the experts
users = readcell(sprintf('%s/users.txt',config.path.users));

% Go through each user
for iuser = 1 : numel(users)
    
    % Load the metadata of the current user
    current_user = users{iuser};
    load(sprintf('%s/dataset/dataset_%s.mat',config.path.users,current_user));
    files = dir(sprintf('%s/trl/%s/*.mat',config.path.users,current_user));
    
    % Go through each file
    for ifile = 1 : numel(files)
        
        % Load the file metadata
        trl = load(sprintf('%s/%s',files(ifile).folder,files(ifile).name));
        
        % Load the raw file
        complete_fname = sprintf('%s/sub-%s/ses-%s/eeg/%s',...
            config.path.curated,trl.subject,trl.session,trl.fileinfo.file);
        
        % Read the header
        header = ft_read_header(complete_fname);
        
        % Gets the data.
        cfg         = [];
        cfg.dataset = complete_fname;
        cfg.header  = header;
        
        wholedata   = my_read_data ( cfg );
        
        
    end
    
end
