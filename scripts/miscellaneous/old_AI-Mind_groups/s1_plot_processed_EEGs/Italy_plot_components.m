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
for ifile = randi(numel(my_dataset)) %1:numel(my_dataset)
    
    % Load the current metadata
    current_metadata = my_dataset(ifile);
    
    fprintf('\n\n %s \n\n',current_metadata.eeg_file)
    
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
    
   % Generates the layout.
    cfg = [];
    cfg.layout = 'eeg1005';
    layout = ft_prepare_layout ( cfg );
    
    % Displays the component data.
    cfg = [];
    cfg.channel         = 1: 15;
    cfg.plotlabels      = 'yes';
    cfg.viewmode        = 'component';
    cfg.layout          = layout;
    cfg.continous       = 'yes';
    cfg.compscale       = 'local';
    cfg.blocksize = 20;
    cfg.linecolor = repmat([217  83  25] / 255,numel(rawdata.label),1);
    
    cfg = ft_databrowser ( cfg, compdata );
    
end

