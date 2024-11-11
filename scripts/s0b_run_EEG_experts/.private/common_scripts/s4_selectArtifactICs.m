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

% Sets the list of component types.
% - Type 0 is always clean component ('Clean component').
% - Type 1 is always EOG component ('EOG component').
% - Type 2 is always EKG component ('EKG component').
config.comptypes      = { 'Clean component' 'EOG component' 'EKG component' 'Noisy component' };

% Sets the visualization configuration parameters.
config.trialfun       = 'restingSegmentation';
config.segment        = 20;
config.overlap        = 10;
config.equal          = false;
config.padding        = 2;
config.addpadd        = false;

config.channel.groups = { 'EEG' };
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
    already_processed = strcmp(metadata.history.previous_steps,'Component revision');
    already_processed = sum(already_processed) > 0;
    if already_processed && ~config.overwrite
        fprintf('  Already processed. Do not overwrite\n\n');
        continue
    end
    
    fprintf ( 1, 'Working on %s.\n', metadata.fileinfo.file );
    
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
        fileinfo            = metadata.fileinfo ( findex );
        artinfo             = metadata.artinfo  ( findex );
        dummy_path_current_file = sprintf('%s/%s',config.path.raw,fileinfo.file);
        
        fprintf ( 1, '    Reading data from disk.\n' );
        
        % Gets the MEG data.
        cfg                   = [];
        cfg.dataset           = dummy_path_current_file;
        cfg.header            = fileinfo.header;
        
        wholedata             = my_read_data ( cfg );
        
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
        fir                   = fir1 ( filtorder, config.filter.band / ( wholedata.fsample / 2 ) );
        wholedata             = myft_filtfilt ( fir, 1, wholedata );
        wholedata             = my_downsample ( wholedata, downrate );
        
        
        % Extracts the overlapping epochs for the artifact revision.
        trialfun              = str2func ( config.trialfun );
        
        fileconfig            = config;
        fileconfig.dataset    = dummy_path_current_file;
        fileconfig.header     = wholedata.hdr;
        fileconfig.begtime    = fileinfo.begtime;
        fileconfig.endtime    = fileinfo.endtime;
        fileconfig.feedback   = 'no';
        
        fileconfig.channel.bad = metadata.chaninfo.bad;
        
        fileconfig.trl        = trialfun ( fileconfig );
        
        trialdata             = ft_redefinetrial ( fileconfig, wholedata );
        clear wholedata;
        
        % Removes the 'cfg' field.
        trialdata             = rmfield ( trialdata, 'cfg' );
        
        trialdata.trial       = cellfun ( @single, trialdata.trial, 'UniformOutput', false );
        
        % Corrects the sample information of the trials.
        trialdata.sampleinfo  = trialdata.sampleinfo + offsets ( findex );
        
        % Resamples and corrects the events definitions.
        event                 = fileinfo.event (:);
        
        if numel ( event )
            [ event.sample ]      = my_deal ( ceil ( [ event.sample ] / downrate ) );
            [ event.sample ]      = my_deal ( [ event.sample ] + offsets ( findex ) );
        end
        
        % Resamples and corrects the artifact definitions.
        artifact              = artinfo.artifact;
        arttypes              = fieldnames ( artifact );
        
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
    event      = cat ( 1, events {:} );
    
    % Merges the artifact definitions.
    artifact   = [];
    arttypes   = fieldnames ( artifacts {1} );
    
    for aindex = 1: numel ( arttypes )
        arttype                       = arttypes { aindex };
        dummy                         = cellfun ( @(x) x.( arttype ).artifact, artifacts, 'UniformOutput', false );
        artifact.( arttype ).artifact = cat ( 1, dummy {:} );
    end
    
    % Removes the undesired events, if requested.
    if config.event < Inf
        event               = event ( ismember ( [ event.value ], config.event ) );
    end
    
    
    fprintf ( 1, '  Combining the components for all channel groups.\n' );
    
    % Selects the SOBI components for the last step.
    compinfo             = metadata.compinfo;
    SOBI                 = compinfo.SOBI;
    
    % Combines the components for all the chanel types.
    channels             = fieldnames ( SOBI );
    channels             = channels ( ismember ( channels, config.channel.groups ) );
    compchans            = cell ( numel ( channels ), 1 );
    compposs             = cell ( numel ( channels ), 1 );
    topolabels           = cell ( numel ( channels ), 1 );
    unmixings            = cell ( numel ( channels ), 1 );
    mixings              = cell ( numel ( channels ), 1 );
    EOGcomps             = cell ( numel ( channels ), 1 );
    EKGcomps             = cell ( numel ( channels ), 1 );
    
    % Goes through each channel type.
    for chindex = 1: numel ( channels )
        
        channel                = channels { chindex };
        
        % Gets the values for the current channel type.
        compchans  { chindex } = chindex * ones ( numel ( SOBI.( channel ).topolabel ), 1 );
        compposs   { chindex } = ( 1: numel ( SOBI.( channel ).topolabel ) )';
        topolabels { chindex } = SOBI.( channel ).topolabel;
        mixings    { chindex } = SOBI.( channel ).mixing;
        unmixings  { chindex } = SOBI.( channel ).unmixing;
        
        % Gets the values for the rejected components.
        if isfield ( SOBI.( channel ), 'type' )
            EOGcomps { chindex } = find ( SOBI.( channel ).type == 1 );
            EOGcomps { chindex } = cat ( 2, chindex * ones ( size ( EOGcomps { chindex } ) ), EOGcomps { chindex } );
            
            EKGcomps { chindex } = find ( SOBI.( channel ).type == 2 );
            EKGcomps { chindex } = cat ( 2, chindex * ones ( size ( EKGcomps { chindex } ) ), EKGcomps { chindex } );
        else
            EOGcomps { chindex } = zeros ( 0, 2 );
            EKGcomps { chindex } = zeros ( 0, 2 );
        end
    end
    
    % Concatenates the matrices.
    compchan            = cat     ( 1, compchans  {:} );
    comppos             = cat     ( 1, compposs   {:} );
    topolabel           = cat     ( 1, topolabels {:} );
    mixing              = blkdiag ( mixings       {:} );
    unmixing            = blkdiag ( unmixings     {:} );
    EOGcomp             = cat     ( 1, EOGcomps   {:} );
    EKGcomp             = cat     ( 1, EKGcomps   {:} );
    
    % Tries to sort the components.
    %%% TO DO
    
    fprintf ( 1, '  Calculating the components.\n' );
    
    % Gets the components for all the channel types at once.
    cfg = [];
    cfg.topolabel       = topolabel;
    cfg.mixing          = mixing;
    cfg.unmixing        = unmixing;
    cfg.updatesens      = 'no';
    cfg.demean          = 'no';
    cfg.doscale         = 'no';
    cfg.feedback        = 'no';
    
    compdata            = my_componentanalysis ( cfg, trialdata );
    
    % Fixes the component names.
    compdata.label      = strrep ( compdata.label, 'component', strcat ( channels ( compchan ), 'C' ) );
    
    % Corrects the topology maps using the original mixing matrix.
    for mychan = 1: numel ( topolabel )
        
        ftchan = strcmp ( compdata.topolabel, topolabel { mychan } );
        
        mytopo = mixing ( mychan, : );
        mytopo ( mytopo == 0 ) = NaN;
        
        compdata.topo ( ftchan, : ) = mytopo;
    end
    
    % Generates the layout.
    layout              = my_prepare_layout ( trialdata );
    compdata.grad       = [];
    compdata.elec       = [];
    
    % Displays the component data.
    cfg = [];
    cfg.channel         = 1: 15;
    cfg.plotlabels      = 'yes';
    cfg.viewmode        = 'component';
    cfg.layout          = layout;
    cfg.continous       = 'yes';
    cfg.compscale       = 'local';
    cfg.colorgroups     = compchan;
    cfg.channelcolormap = [   0 114 189;   0 114 189; 162  20  47; 217  83  25; 126  47 142 ] / 255;
    cfg.artfctdef       = artifact;
    cfg.ploteventlabels = 'colorvalue';
    cfg.event           = event;
    
    my_databrowser ( cfg, compdata );
    
    figures             = gcf;
    
    
    compdata.topo  ( isnan ( compdata.topo ) ) = 0;
    
    % Calculates the position of the current EOG and EKG components.
    EOGcomp             = find ( ismember ( [ compchan comppos ], EOGcomp, 'rows' ) );
    EKGcomp             = find ( ismember ( [ compchan comppos ], EKGcomp, 'rows' ) );
    
    % Iterates until the artifact components are correctly removed.
    while true
        
        % Asks for the EOG and EKG components to remove.
        question            = { 'EOG components:' 'EKG components:' };
        default             = { num2str( EOGcomp' ) num2str( EKGcomp' ) };
        answer              = mydlg_inputdlg ( question, 'Question', 1, default );
        
        if isempty ( answer ), break, end
        
        EOGcomp             = str2num ( answer {1} ); %#ok<ST2NM>
        EKGcomp             = str2num ( answer {2} ); %#ok<ST2NM>
        
        EOGcomp             = EOGcomp (:);
        EKGcomp             = EKGcomp (:);
        
        % Removes the invalid components.
        EOGcomp ( EOGcomp < 1 ) = [];
        EOGcomp ( EOGcomp > numel ( compdata.label ) ) = [];
        EKGcomp ( EKGcomp < 1 ) = [];
        EKGcomp ( EKGcomp > numel ( compdata.label ) ) = [];
        
        EOGcomp             = EOGcomp (:);
        EKGcomp             = EKGcomp (:);
        
        % If a IC appears in both groups rises an error.
        if numel ( intersect ( EOGcomp, EKGcomp ) )
            
            uiwait ( errordlg ( 'Some of the selected components are labelled as both EOG and EKG.', 'Error' ) );
            continue
        end
        
        drawnow
        
        % Reconstructs the data without the selected components.
        cfg                 = [];
        cfg.component       = cat ( 1, EOGcomp, EKGcomp );
        cfg.demean          = 'no';
        
        cleandata           = my_rejectcomponent ( cfg, compdata, trialdata );
        
        
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
        cfg                 = [];
        cfg.channel         = 1: 32;
        cfg.physio          = physio;
        cfg.gradscale       = 0.05;
        cfg.eegscale        = 1e-13;
        cfg.eogscale        = 1e-13;
        cfg.ecgscale        = 1e-13;
        cfg.ylim            = [ -1 1 ] * 2.5e-12;
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
        
        
        my_databrowser ( cfg, trialdata );
        figures (2)         = gcf;
        my_databrowser ( cfg, cleandata );
        figures (3)         = gcf;
        drawnow
        
        answer              = mydlg_questdlg ( 'Are the EOG and EKG artifacts correctly removed?', 'Question', 'Yes', 'No', 'Cancel', 'Yes' );
        delete ( figures ( 2: 3 ) )
        
        if ~strcmp ( answer, 'No' ), break, end
    end
    
    % Deletes the data browser figures.
    delete ( intersect ( figures, findall (0) ) )
    drawnow
    
    if strcmp  ( answer, 'Cancel' ), return, end
    if isempty ( answer ),           return, end
    
    
    % Stores each component to remove in its channel type.
    for chindex = 1: numel ( channels )
        
        channel  = channels { chindex };
        
        % Gets the EOG and EKG components for the current channel type.
        cEOGcomp = comppos ( EOGcomp ( compchan ( EOGcomp ) == chindex ) );
        cEKGcomp = comppos ( EKGcomp ( compchan ( EKGcomp ) == chindex ) );
        
        % Labels the EOG and EKG components.
        comptype = zeros ( size ( SOBI.( channel ).topolabel ) );
        comptype ( cEOGcomp ) = 1;
        comptype ( cEKGcomp ) = 2;
        
        SOBI.( channel ).type = comptype;
    end
    
    
    % Updates the current step structure.
    compinfo.step        = 'Component revision';
    compinfo.date        = datestr ( now );
    compinfo.config      = config;
    compinfo.types       = config.comptypes;
    compinfo.SOBI        = SOBI;
    metadata.compinfo    = compinfo;
    
    % history
    metadata.history.last_step = 'Component revision';
    metadata.history.previous_steps(end+1) = {metadata.history.last_step};
    metadata.history.date = datestr ( now );
    
    % Saves the output data.
    save ( '-v6', metadata_file, '-struct', 'metadata' )
end
