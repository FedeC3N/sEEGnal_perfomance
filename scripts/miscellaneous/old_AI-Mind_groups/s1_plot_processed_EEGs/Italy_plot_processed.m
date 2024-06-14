clear
clc
close all

% Paths
path = [];
path.metadata = '../../data/metadata/unified_format';
path.dataset = '../../data/metadata/dataset';

% Add fieldtrip
addpath('../../../../SharedFunctions/fieldtrip-20220626');
addpath('../../../../SharedFunctions/functions');
ft_defaults

% Files (based on dataset)
load(sprintf('%s/dataset_all.mat',path.dataset));

% Keep only Fede_processed files
files_included_mask = strcmp({my_dataset.processed_by},'Italy_processed');
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
    % Italy did not exclude badchannels for the ICA estimation so we can
    % not exclude badchannels either.
    badchannels = metadata.badchannels.badchannels_labels;
    included_channels_mask = exclude_physio_channel;
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
    
    % Define the artifacts as visuals
    visual = [];
    if numel(metadata.annotations.start_sample) == 0
        visual.artifact = zeros(1,2);
    else
        visual.artifact = cat(2,metadata.annotations.start_sample,...
            metadata.annotations.end_sample);
        % Have to apply the downsampling
        downsampling_rate = new_fs/original_fs;
        visual.artifact = downsampling_rate*visual.artifact;
        visual.artifact = ceil(visual.artifact);
        
    end
    
    % Databrowser
    cfg = [];
    cfg.channel = 1:32;
    cfg.continous       = 'yes';
    cfg.blocksize = 20;
%     cfg.ylim = [-50000 50000];
    cfg.viewmode = 'vertical';
    cfg.artfctdef.visual.artifact     = visual.artifact;
    cfg.linecolor = repmat([217  83  25] / 255,numel(rawdata.label),1);
    cfg.artifactalpha = 0.5;
    cfg = ft_databrowser ( cfg, cleandata );
    
end

