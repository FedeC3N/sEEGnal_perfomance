clc
clear
close all

% Adds the functions folders to the path.
addpath ( '../../SharedFunctions/fieldtrip-20220626/');
addpath ( '../../SharedFunctions/functions/');
addpath ( '../../SharedFunctions/mne_silent/');
addpath ( '../../SharedFunctions/functions_eep/');
ft_defaults

% Select user
fid = fopen('./user.txt','r');
user = fgetl(fid);
fclose(fid);

% Paths
config.path.raw  = '../../../data/raw' ;
config.path.trl  = ['../../../meta/trl/' user];
config.path.dataset_info = ['../../../meta/dataset/dataset_' user '.mat'];

% Action when the task have already been processed.
config.overwrite      = false;

% Sets the segmentation parameters.
config.trialfun       = 'restingSegmentation';
config.segment        = 4;
config.padding        = 2;
config.addpadd        = false;

% Artifacts to exclude of the IC analysis.
config.artifact       = { 'visual' 'jump' 'muscle' };

% Sets the IC analisys parameters.
config.channel.groups = { 'EEG' };
config.channel.ignore = {};

% Sets the filter band.
config.filter.band    = [ 2 45 ];

% Lists the files.
dataset    = struct2array ( load ( config.path.dataset_info ) );

% Goes through each subject and task.
for sindex =  1: numel ( dataset )
    
    % Read the metadata
    current_file = dataset(sindex).file;
    current_subject = dataset(sindex).subject;
    current_session = dataset(sindex).session;
    current_task = dataset(sindex).task;
    fprintf('Current file sub-%s_ses-%s_task-%s\n', current_subject,current_session,current_task);
    
    % Check if exist metadata file (mandatory)
    metadata_file   =  sprintf ( '%s/sub-%s_ses-%s_task-%s_eeg.mat', config.path.trl, ...
        current_subject, current_session, current_task );
    if ~exist(metadata_file)
        fprintf('Metadata file not created. Skip \n\n');
    end
    
    % Load metadata
    metadata = load (metadata_file);
    dummy_path_current_file = sprintf('%s/%s',config.path.raw,current_file);
    
    % Check if overwrite
    already_processed = strcmp(metadata.history.previous_steps,'Component extraction');
    already_processed = sum(already_processed) > 0;
    if already_processed && ~config.overwrite
        fprintf('  Already processed. Do not overwrite\n\n');
        continue
    end
    
    fprintf ( 1, 'Working on %s.\n', metadata.fileinfo.file );
    
    % Reserves memory for the data segments.
    datas             = cell ( numel ( metadata.fileinfo ), 1 );
    headers           = cell ( numel ( metadata.fileinfo ), 1 );
    
    % Calculates the optimal downsampling rate for the frequency band.
    downrate          = floor ( metadata.fileinfo (1).header.Fs / ( 2 * config.filter.band ( end ) ) );
    
    % Calculates the optimal filter order from the desired padding.
    filtorder         = floor ( metadata.fileinfo (1).header.Fs * config.padding );
    
    
    % Goes through all the data files.
    for findex = 1: numel ( metadata.fileinfo )
        
        % Gets the current epochs and artifact definitions.
        fileinfo            = metadata.fileinfo  ( findex );
        artinfo             = metadata.artinfo ( findex );
        dummy_path_current_file = sprintf('%s/%s',config.path.raw,fileinfo.file);

        fprintf ( 1, '    Reading data from disk.\n' );
        
        % Gets the MEG data.
        cfg                 = [];
        cfg.dataset         = dummy_path_current_file;
        cfg.header          = fileinfo.header;
        
        wholedata           = my_read_data ( cfg );
        
        % Selects the channels.
        cfg                  = [];
        badchannels          = strcat('-',metadata.chaninfo.bad);
        desired_channel      = cat(2,{'eeg'},badchannels);
        cfg.channel          = desired_channel;
        cfg.precision        = 'single';
        cfg.feedback         = 'no';
        
        wholedata            = ft_preprocessing ( cfg, wholedata );
        
        
        fprintf ( 1, '    Filtering the data in the band %0.0f - %0.0f Hz.\n', config.filter.band );
        
        % Filters and downsamples the data.
        fir                 = fir1 ( filtorder, config.filter.band / ( wholedata.fsample / 2 ) );
        wholedata           = myft_filtfilt ( fir, 1, wholedata );
        wholedata           = my_downsample ( wholedata, downrate );
        
        
        % Resamples and corrects the artifact definitions.
        artifact            = artinfo.artifact;
        arttypes            = fieldnames ( artifact );
        
        for aindex = 1: numel ( arttypes )
            arttype                       = arttypes { aindex };
            artifact.( arttype ).artifact = ceil ( artifact.( arttype ).artifact / downrate );
        end
        
        % Gets the artifact free epochs.
        trialfun            = str2func ( config.trialfun );
        
        fileconfig          = config;
        fileconfig.dataset  = dummy_path_current_file;
        fileconfig.header   = wholedata.hdr;
        fileconfig.begtime  = fileinfo.begtime;
        fileconfig.endtime  = fileinfo.endtime;
        fileconfig.feedback = 'no';
        
        fileconfig.artifact = artifact;
        fileconfig.artifact = rmfield ( fileconfig.artifact, setdiff ( fieldnames ( artinfo.artifact ), config.artifact ) );
        
        fileconfig.channel.bad = metadata.chaninfo.bad;
        
        fileconfig.trl      = trialfun ( fileconfig );
        
        trialdata           = ft_redefinetrial ( fileconfig, wholedata );
        trialdata.trial     = cellfun ( @single, trialdata.trial, 'UniformOutput', false );
        
        
        % Stores the epoch data.
        datas   { findex }  = cat ( 3, trialdata.trial {:} );
        headers { findex }  = trialdata.hdr;
        
        clear trialdata
    end
    
    % Converts the data cell to a matrix.
    datas             = cat ( 3, datas {:} );
    
    if isempty ( datas )
        fprintf ( 1, '  Ignoring subject ''%s'' (no clean data found in any file).\n', metadata.subject );
        continue
    end
    
    % Limits the number of trials to a maximum of 200.
    if size ( datas, 3 ) > 200
        sample = sort ( randsample ( size ( datas, 3 ), 200 ) );
        datas  = datas ( :, :, sample );
    end
    
    
    % Recovers the SOBI information.
    if isfield ( metadata, 'compinfo' )
        compinfo          = metadata.compinfo;
        
        SOBI              = compinfo.SOBI;
        
    % If no SOBI information initializes the information structure.
    else
        compinfo          = [];
        compinfo.step     = [];
        compinfo.date     = [];
        compinfo.config   = [];
        compinfo.types    = [];
        compinfo.SOBI     = [];
        compinfo.history  = {};
        
        SOBI              = [];
    end
    
    % Goes through each channel group.
    for chindex = 1: numel ( config.channel.groups )
        
        channel  = config.channel.groups { chindex };
        
        % Gets the labels for this channel group.
        label    = ft_channelselection ( channel, headers {1}.label );
        
        % Ignores the selected channels.
        label    = setdiff ( label, config.channel.ignore );
        label    = setdiff ( label, metadata.chaninfo.bad );
        
        if isempty ( label )
            fprintf ( 1, '  Ignoring channel group ''%s'' (no data).\n', channel );
            continue
        end
        
        fprintf ( 1, '  Working with channel group ''%s''.\n', channel );
        
        % Gets the channels and labels in the right order.
        chanindx = ismember ( headers {1}.label, label );
        label    = headers {1}.label ( chanindx );
        chandata = datas ( chanindx, :, : );
        
        fprintf ( 1, '    Extracting the SOBI components.\n' );
        fprintf ( 1, '      ' );
        
        % Gets the independen components for the current channel group.
        mixing   = my_sobi ( chandata );
        
        % Adds the channel group to the output.
        SOBI.( config.channel.groups { chindex } )           = [];
        SOBI.( config.channel.groups { chindex } ).topolabel = label;
        SOBI.( config.channel.groups { chindex } ).mixing    = mixing;
        SOBI.( config.channel.groups { chindex } ).unmixing  = pinv ( mixing );
        SOBI.( config.channel.groups { chindex } ).type      = zeros ( size ( label ) );
    end
    
    % Updates the current step structure.
    compinfo.config   = config;
    compinfo.types    = { 'Clean component' };
    compinfo.SOBI     = SOBI;
    
    % Adds the SOBI information to the task information.
    metadata.compinfo = compinfo;
    
    % history
    metadata.history.last_step = 'Component extraction';
    metadata.history.previous_steps(end+1) = {metadata.history.last_step};
    metadata.history.date = datestr ( now );
    
    % Saves the intependent components.
    save ( '-v6', metadata_file, '-struct', 'metadata' )
    
end
