clear
clc
close all

% Paths
path = [];
path.metadata = '../../data/metadata/unified_format';
path.dataset = '../../data/metadata/dataset';

% Add fieldtrip
addpath('../../../../SharedFunctions/fieldtrip-20220626');
addpath('../../../../SharedFunctions/functions/');
ft_defaults

% Files (based on dataset)
load(sprintf('%s/dataset_all.mat',path.dataset));

% Keep only Fede_processed files
files_included_mask = strcmp({my_dataset.processed_by},'Fede_processed');
my_dataset = my_dataset(files_included_mask);

% For each file
for ifile = 1:numel(my_dataset)
    
    % Load the current metadata
    current_metadata = my_dataset(ifile);
  
    % Load the metadata and extract the information
    metadata = load(sprintf('%s/%s',current_metadata.metadata_path,...
        current_metadata.metadata_file));
    
    % Load the EEG recording
    cfg = [];
    dummy_path_current_file = sprintf('%s/%s',current_metadata.eeg_path,...
        current_metadata.eeg_file);
    cfg.dataset = dummy_path_current_file;
    cfg.continuous = 'yes';
    % - filtered 2-45 Hz
    cfg.bpfilter = 'yes';
    cfg.bpfreq = [2 45];
    % - Exclude channels
    all_channels = metadata.badchannels.channel_labels;
    physio_channels = {'VEOGL', 'CLAV','EMG1','EMG2'};
    exclude_physio_channel = ~ismember(all_channels,physio_channels);
    % - Exclude badchannels
    badchannels = metadata.badchannels.badchannels_labels;
    if isempty(badchannels)
        included_channels_mask = exclude_physio_channel;
    else
        exclude_badchannels =  ~ismember(all_channels,badchannels);
    included_channels_mask = exclude_physio_channel & exclude_badchannels;
    end    
    % Channel selection
    cfg.channel = all_channels(included_channels_mask);
    % - Re-reference to commmon average
    cfg.reref = 'yes';
    cfg.refchannel = 'all';
    rawdata = ft_preprocessing (cfg);   
    
    % Downsample
    original_fs = rawdata.fsample;
    new_fs = 100;
    cfg = [];
    cfg.resamplefs = new_fs;
    rawdata = ft_resampledata(cfg, rawdata);    
    
    % Obtain the component data
    % Create the compinformation
    mixing = metadata.ICs.mixing_matrix;
    unmixing = metadata.ICs.unmixing_matrix;
    topolabel = all_channels(included_channels_mask);
    excluded_components = find(metadata.ICs.IC_included_mask == 0);
    
    % Gets the components for all the channel types at once.
    cfg = [];
    cfg.topolabel       = topolabel;
    cfg.mixing          = mixing;
    cfg.unmixing        = unmixing;
    cfg.updatesens      = 'no';
    cfg.demean          = 'no';
    cfg.doscale         = 'no';
    cfg.feedback        = 'no';
    
    compdata            = my_componentanalysis ( cfg, rawdata );
    
    % Change the labels
    compdata.label      = strrep ( compdata.label, 'component', 'IC');
    
    % Corrects the topology maps using the original mixing matrix.
    for mychan = 1: numel ( topolabel )
        
        ftchan = strcmp ( compdata.topolabel, topolabel { mychan } );
        
        mytopo = mixing ( mychan, : );
        mytopo ( mytopo == 0 ) = NaN;
        
        compdata.topo ( ftchan, : ) = mytopo;
    end
    
    % Reject artifacts components
    cfg                 = [];
    cfg.component       = excluded_components;
    cfg.demean          = 'no';
    cleandata           = my_rejectcomponent ( cfg, compdata, rawdata );
    
    % Define the artifacts
    eog = [];
    eog_index = strcmp(metadata.annotations.artifact_type,'EOG');
    if sum(eog_index) == 0
        eog.artifact = zeros(0,2);
    else
        start_sample = metadata.annotations.start_sample(eog_index);
        end_sample = metadata.annotations.end_sample(eog_index);
        eog.artifact = cat(2,start_sample,end_sample);
        % Have to apply the downsampling
        downsampling_rate = new_fs/original_fs;
        eog.artifact = downsampling_rate*eog.artifact;
    end
    muscle = [];
    muscle_index = strcmp(metadata.annotations.artifact_type,'Muscle');
    if sum(muscle_index) == 0
        muscle.artifact = zeros(0,2);
    else
        start_sample = metadata.annotations.start_sample(muscle_index);
        end_sample = metadata.annotations.end_sample(muscle_index);
        muscle.artifact = cat(2,start_sample,end_sample);
        % Have to apply the downsampling
        downsampling_rate = new_fs/original_fs;
        muscle.artifact = downsampling_rate*muscle.artifact;
    end
    jump = [];
    jump_index = strcmp(metadata.annotations.artifact_type,'Jump');
    if sum(jump_index) == 0
        jump.artifact = zeros(0,2);
    else
        start_sample = metadata.annotations.start_sample(jump_index);
        end_sample = metadata.annotations.end_sample(jump_index);
        jump.artifact = cat(2,start_sample,end_sample);
        % Have to apply the downsampling
        downsampling_rate = new_fs/original_fs;
        jump.artifact = downsampling_rate*jump.artifact;
    end
    visual = [];
    visual_index = strcmp(metadata.annotations.artifact_type,'Visual');
    if sum(visual_index) == 0
        visual.artifact = zeros(0,2);
    else
        start_sample = metadata.annotations.start_sample(visual_index);
        end_sample = metadata.annotations.end_sample(visual_index);
        visual.artifact = cat(2,start_sample,end_sample);
        % Have to apply the downsampling
        downsampling_rate = new_fs/original_fs;
        visual.artifact = downsampling_rate*visual.artifact;
    end
    
    % Databrowser
    cfg = [];
    cfg.channel = 1:32;
    cfg.continous       = 'yes';
    cfg.blocksize = 20;
    cfg.ylim = [-30 30];
    cfg.viewmode = 'vertical';
    cfg.channelclamped = {'EOG'}; 
    cfg.artfctdef.eog.artifact     = eog.artifact;
    cfg.artfctdef.muscle.artifact     = muscle.artifact;
    cfg.artfctdef.jump.artifact     = jump.artifact;
    cfg.artfctdef.visual.artifact     = visual.artifact;
    cfg.linecolor = repmat([217  83  25] / 255,numel(rawdata.label),1);
    cfg.artifactalpha = 0.5;
    cfg = ft_databrowser ( cfg, cleandata );
    
    % Construct the metadata
    eog_label = cell(size(cfg.artfctdef.eog.artifact,1),1);
    eog_label(:) = {'EOG'};
    muscle_label = cell(size(cfg.artfctdef.muscle.artifact,1),1);
    muscle_label(:) = {'Muscle'};
    jump_label = cell(size(cfg.artfctdef.jump.artifact,1),1);
    jump_label(:) = {'Jump'};
    visual_label = cell(size(cfg.artfctdef.visual.artifact,1),1);
    visual_label(:) = {'Visual'};
    artifact_types = cat(1,eog_label,muscle_label,jump_label,visual_label);
    start_sample = cat(1,cfg.artfctdef.eog.artifact(:,1),cfg.artfctdef.muscle.artifact(:,1),...
        cfg.artfctdef.jump.artifact(:,1),cfg.artfctdef.visual.artifact(:,1));
    start_sample = start_sample * original_fs/new_fs; %Upsample
    end_sample = cat(1,cfg.artfctdef.eog.artifact(:,2),cfg.artfctdef.muscle.artifact(:,2),...
        cfg.artfctdef.jump.artifact(:,2),cfg.artfctdef.visual.artifact(:,2));
    end_sample = end_sample * original_fs/new_fs; %Upsample
    
    % Save metadata
    metadata.annotations.artifact_type = artifact_types;
    metadata.annotations.start_sample = start_sample;
    metadata.annotations.end_sample = end_sample;
    
    save(sprintf('%s/%s',current_metadata.metadata_path,...
        current_metadata.metadata_file),'-struct','metadata')
    
    
end

