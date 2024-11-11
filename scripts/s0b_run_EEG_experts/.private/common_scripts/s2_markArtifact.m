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

% Sets the visualization configuration parameters.
config.trialfun       = 'restingSegmentation';
config.segment        = 20;
config.overlap        = 10;
config.equal          = false;
config.padding        = 2;
config.addpadd        = false;

config.channel.data   = {'EEG' };
config.channel.ignore = {};

% Sets the filter band.
config.filter.band    = [ 2 45 ];

% Sets the events to show (Inf for all).
config.event          = Inf;

% Lists the files.
dataset    = struct2array ( load ( config.path.dataset_info ) );

% Goes through each subject and task.
for sindex = 1: numel ( dataset )
    
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
    already_processed = strcmp(metadata.history.previous_steps,'Select artifacts');
    already_processed = sum(already_processed) > 0;
    if already_processed && ~config.overwrite
        fprintf('  Already processed. Do not overwrite\n\n');
        continue
    end
    
    fprintf ( 1, 'Working on %s.\n', current_file );
    
    % Gets the files length as a vector.
    headers             = [ metadata.fileinfo.header ];
    samples             = [ headers.nSamples ];
    offsets             = cat ( 2, 0, cumsum ( samples) );
    
    % Calculates the optimal downsampling rate for the frequency band.
    downrate            = floor ( headers (1).Fs / ( 2 * config.filter.band ( end ) ) );
    offsets             = ceil  ( offsets / downrate );
    
    % Calculates the optimal filter order from the desired padding.
    filtorder           = floor ( headers (1).Fs * config.padding );
    
    
    % Reserves memory for the trialdata and the artifact definitions.
    trialdatas          = cell ( numel ( metadata.fileinfo ), 1 );
    events              = cell ( numel ( metadata.fileinfo ), 1 );
    artifacts           = cell ( numel ( metadata.fileinfo ), 1 );
    
    % Goes through all the data files.
    for findex = 1: numel ( metadata.fileinfo )
        
        % Gets the current epochs and artifact definitions.
        fileinfo             = metadata.fileinfo ( findex );
        artinfo              = metadata.artinfo  ( findex );
        dummy_path_current_file = sprintf('%s/%s',config.path.raw,fileinfo.file);
        
        fprintf ( 1, '    Reading data from disk.\n' );
        
        % Gets the MEG data.
        cfg                  = [];
        cfg.dataset          = dummy_path_current_file;
        cfg.header           = fileinfo.header;
        
        wholedata            = my_read_data ( cfg );
        
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
        fir                  = fir1 ( filtorder, config.filter.band / ( wholedata.fsample / 2 ) );
        wholedata            = myft_filtfilt ( fir, 1, wholedata );
        wholedata            = my_downsample ( wholedata, downrate );
        
        
        % Extracts the overlapping epochs for the artifact revision.
        trialfun             = str2func ( config.trialfun );
        
        fileconfig           = config;
        fileconfig.dataset   = dummy_path_current_file;
        fileconfig.header    = wholedata.hdr;
        fileconfig.begtime   = fileinfo.begtime;
        fileconfig.endtime   = fileinfo.endtime;
        fileconfig.feedback  = 'no';
        
        fileconfig.channel.bad = metadata.chaninfo.bad;
        
        fileconfig.trl       = trialfun ( fileconfig );
        
        trialdata            = ft_redefinetrial ( fileconfig, wholedata );
        clear wholedata;
        
        % Removes the 'cfg' field.
        trialdata            = rmfield ( trialdata, 'cfg' );
        
        trialdata.trial      = cellfun ( @single, trialdata.trial, 'UniformOutput', false );
        
        % Corrects the sample information of the trials.
        trialdata.sampleinfo = trialdata.sampleinfo + offsets ( findex );
        
        % Resamples and corrects the event definitions.
        event                = fileinfo.event (:);
        
        if numel ( event )
            [ event.sample ]      = my_deal ( ceil ( [ event.sample ] / downrate ) );
            [ event.sample ]      = my_deal ( [ event.sample ] + offsets ( findex ) );
        end
        
        % Resamples and corrects the artifact definitions.
        artifact             = artinfo.artifact;
        arttypes             = fieldnames ( artifact );
        
        for aindex = 1: numel ( arttypes )
            arttype                       = arttypes { aindex };
            artifact.( arttype ).artifact = ceil ( artifact.( arttype ).artifact / downrate );
            artifact.( arttype ).artifact = artifact.( arttype ).artifact + offsets ( findex );
        end
        
        
        % Stores the epoch data.
        trialdatas  { findex } = trialdata;
        events      { findex } = event;
        artifacts   { findex } = artifact;
    end
    
    if all ( cellfun ( @isempty, trialdatas ) )
        fprintf ( 1, '  No data. Ignoring.\n' );
        continue
    end
    
    fprintf ( 1, '  Merging all the data files.\n' );
    
    % Merges the epoch data.
    cfg                 = [];
    
    trialdata           = myft_appenddata ( cfg, trialdatas {:} );
    clear trialdatas
    
    % Merges the event information.
    event               = cat ( 1, events {:} );
    
    % Merges the artifact definitions.
    artifact            = [];
    arttypes            = fieldnames ( artifacts {1} );
    
    for aindex = 1: numel ( arttypes )
        arttype                       = arttypes { aindex };
        dummy                         = cellfun ( @(x) x.( arttype ).artifact, artifacts, 'UniformOutput', false );
        artifact.( arttype ).artifact = cat ( 1, dummy {:} );
    end
    
    % Removes the undesired events, if requested.
    if config.event < Inf
        event               = event ( ismember ( [ event.value ], config.event ) );
    end
    
    
    % Gets the list of non-MEG channels (EOG, EKG, etc.).
    MEG                 = find ( ismember ( trialdata.label', ft_channelselection ( { 'MEG' },             trialdata.label ) ) );
    MEGMAG              = find ( ismember ( trialdata.label', ft_channelselection ( { 'MEGMAG' },          trialdata.label ) ) );
    MEGGRAD             = find ( ismember ( trialdata.label', ft_channelselection ( { 'MEGGRAD' },         trialdata.label ) ) );
    EEG                 = find ( ismember ( trialdata.label', ft_channelselection ( { 'EEG' },             trialdata.label ) ) );
    physio              = find ( ismember ( trialdata.label', ft_channelselection ( { 'EOG' 'ECG' 'EMG' }, trialdata.label ) ) );
    
    % Sets the colors according to the channel type.
    chancol             = zeros ( numel ( trialdata.label ), 1 );
    chancol ( MEG )     = 1;
    chancol ( MEGMAG )  = 2;
    chancol ( MEGGRAD ) = 3;
    chancol ( EEG )     = 4;
    chancol ( physio )  = 5;
    
    % Defines the EEG projector.
    projector           = [];
    projector.kind      = 10; % EEG average reference.
    projector.active    = true;
    projector.desc      = 'EEG average reference';
    projector.data      = [];
    projector.data.nrow      = 1;
    projector.data.ncol      = sum ( EEG );
    projector.data.row_names = [];
    projector.data.col_names = trialdata.label ( EEG );
    projector.data.data      = ones ( 1, sum ( EEG ) );
    
    % Displays the MEG data to append or remove artifacts.
    cfg = [];
    cfg.channel         = 1: 32;
    cfg.physio          = physio;
    cfg.eegscale        = 1e-13;
    cfg.eogscale        = 1e-13;
    cfg.ecgscale        = 1e-13;
    cfg.ylim            = [ -1 1 ] * 4e-12;
    cfg.plotlabels      = 'yes';
    cfg.viewmode        = 'vertical';
    cfg.continous       = 'yes';
    cfg.colorgroups     = chancol;
    cfg.channelcolormap = [   0 114 189;   0 114 189; 162  20  47; 217  83  25; 126  47 142 ] / 255;
    cfg.artfctdef       = artifact;
    cfg.ploteventlabels = 'colorvalue';
    cfg.event           = event;
    
    cfg.badchan         = fileconfig.channel.bad;
    cfg.projinfo        = projector;
    cfg.applyprojector  = true;
    
    cfg                 = my_databrowser ( cfg, trialdata );
    drawnow
    
    % Goes through all the data files.
    for findex = 1: numel ( metadata.artinfo )
        
        % Uses the full list of artifacts.
        artifact = cfg.artfctdef;
        
        % Gets the previous step artifact definition.
        artinfo         = metadata.artinfo ( findex );
        
        % Gets the number of previous and current samples.
        psamples        = sum ( [ headers( 1: findex - 1 ).nSamples ] );
        csamples        = headers( findex ).nSamples;
        
        % Corrects the artifact definitions.
        arttypes = fieldnames ( artifact );
        for aindex = 1: numel ( arttypes )
            
            % Gets the list of artifacts of the current type.
            arttype          = arttypes { aindex };
            artifact         = cfg.artfctdef.( arttype ).artifact;
            
            % Resamples and corrects the artifact definitions.
            artifact         = artifact * downrate;
            artifact         = artifact - psamples;
            
            % Removes the artifacts of previous files.
            artifact ( artifact ( :, 2 ) < 1,        : ) = [];
            
            % Removes the artifacts of future files.
            artifact ( artifact ( :, 1 ) > csamples, : ) = [];
            
            % Corrects the first and last artifact, if needed.
            artifact = max ( artifact, 1 );
            artifact = min ( artifact, csamples );
            
            % Stores the current artifact definition.
            artinfo.artifact.( arttype ).artifact = artifact;
        end
        
        % Updates the artifact info.
        metadata.artinfo ( findex ) = artinfo;
        
    end
    
    % history
    metadata.history.last_step = 'Select artifacts';
    metadata.history.previous_steps(end+1) = {metadata.history.last_step};
    metadata.history.date = datestr ( now );
    
    % Saves the output data.
    save ( '-v6', metadata_file, '-struct', 'metadata' )
    
end
