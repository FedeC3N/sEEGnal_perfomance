%{

Estimate the variance explained by each individual independent component

@author: Fede

%}

clc
clear
close all

% Adds the functions folders to the path.
addpath ( '../SharedFunctions/fieldtrip-20220626/');
addpath ( '../SharedFunctions/functions/');
addpath ( '../SharedFunctions/mne_silent/');
addpath ( '../SharedFunctions/functions_eep/');
ft_defaults

% Select user
users = {'fede','isa','luis','maria'};

% Paths
config.path.raw  = '../../databases/AI_Mind_database';

% Sets the visualization configuration parameters.
config.trialfun       = 'restingSegmentation';
config.segment        = 20;
config.overlap        = 10;
config.equal          = false;
config.padding        = 2;
config.addpadd        = false;

config.channel.data   = { 'EEG' };
config.channel.ignore = {};

% Sets the filter band.
config.filter.band    = [ 2 45 ];

% Sets the events to show (Inf for all).
original_variance_all = nan(20,4);
reconstructed_variance_all = nan(20,4);
ratio_all = nan(20,4);
for iuser = 1 : numel(users)

    current_user = users{iuser};

    % Lists the files.
    config.path.derivatives = sprintf('../../databases/AI_Mind_database/derivatives/%s',current_user);
    dataset_file = sprintf('%s/dataset_%s.mat',config.path.derivatives,current_user);
    dataset    = struct2array ( load ( dataset_file ) );

    % Goes through each subject and task.
    for sindex = 1: numel ( dataset )

        % Read the metadata
        current_folder = dataset(sindex).folder;
        current_file = dataset(sindex).file;
        current_subject = dataset(sindex).subject;
        current_session = dataset(sindex).session;
        current_task = dataset(sindex).task;
        fprintf('Current file sub-%s_ses-%s_task-%s\n', current_subject,current_session,current_task);

        % Check if exist metadata file (mandatory)
        current_path_trl = sprintf('%s/sub-%s/ses-%s/eeg',...
            config.path.derivatives,current_subject,current_session);
        metadata_file   =  sprintf ( '%s/sub-%s_ses-%s_task-%s_desc-trl_eeg.mat', current_path_trl, ...
            current_subject, current_session, current_task );
        if ~exist(metadata_file)
            fprintf('Metadata file not created. Skip \n\n');
        end

        % Load metadata
        metadata = load (metadata_file);
        dummy_path_current_file = sprintf('%s/%s',current_folder,current_file);

        fprintf ( 1, 'Working on %s.\n', metadata.fileinfo.file );

        fprintf ( 1, '    Reading data from disk.\n' );

        % Gets the MEG data.
        cfg                  = [];
        cfg.dataset          = dummy_path_current_file;
        cfg.header           = metadata.fileinfo.header;

        wholedata            = my_read_data ( cfg );

        % Selects the channels.
        cfg                  = [];
        badchannels          = strcat('-',metadata.chaninfo.bad);
        desired_channel      = cat(2,{'eeg'},badchannels);
        cfg.channel          = desired_channel;
        cfg.precision        = 'single';
        cfg.feedback         = 'no';

        wholedata            = ft_preprocessing ( cfg, wholedata );

        % Filters and downsamples the data.
        downrate            = floor ( metadata.fileinfo.header.Fs / ( 2 * config.filter.band ( end ) ) );
        filtorder           = floor ( metadata.fileinfo.header.Fs * config.padding );
        fir                   = fir1 ( filtorder, config.filter.band / ( wholedata.fsample / 2 ) );
        wholedata             = myft_filtfilt ( fir, 1, wholedata );
        wholedata             = my_downsample ( wholedata, downrate );

        % Estimate the total variance of the original data
        wholedata_matrix = cat(2,wholedata.trial{:});
        original_variance = var(wholedata_matrix,[],2);
        original_variance = sum(original_variance);

        % Get the mixing and unmixing
        SOBI                = metadata.compinfo.SOBI;
        mixing              = SOBI.EEG.mixing;
        unmixing            = SOBI.EEG.unmixing;


        % Removes the artifact components.
        EOGcomp             = find(SOBI.EEG.type == 1);
        EKGcomp             = find(SOBI.EEG.type == 2);
        components_to_reject = cat ( 1, EOGcomp, EKGcomp );

        % Estimate the new variance
        sources_matrix = unmixing * wholedata_matrix;
        sources_matrix(components_to_reject,:) = [];
        mixing(:,components_to_reject) = [];
        reconstructed_matrix = mixing * sources_matrix;
        reconstructed_variance = var(reconstructed_matrix,[],2);
        reconstructed_variance = sum(reconstructed_variance);

        original_variance_all(sindex,iuser) = original_variance;
        reconstructed_variance_all(sindex,iuser) = reconstructed_variance;
        ratio_all(sindex,iuser) = reconstructed_variance/original_variance;

        % Displays the MEG data to append or remove artifacts.
        % cfg                 = [];
        % cfg.channel         = 1: 30;
        % cfg.gradscale       = 0.05;
        % cfg.eegscale        = 1e-13;
        % cfg.eogscale        = 1e-13;
        % cfg.ecgscale        = 1e-13;
        % cfg.ylim            = [ -1 1 ] * 4e-12;
        % cfg.plotlabels      = 'yes';
        % cfg.viewmode        = 'vertical';
        % cfg.continous       = 'yes';
        % cfg.channelcolormap = [   0 114 189;   0 114 189; 162  20  47; 217  83  25; 126  47 142 ] / 255;
        % cfg.ploteventlabels = 'colorvalue';
        %
        % cfg                 = my_databrowser ( cfg, reconstruted );
        % drawnow

        % Saves the output data.
        % metadata.compinfo.SOBI.EEG.explained_variance = explained_variance;
        % save ( '-v6', metadata_file, '-struct', 'metadata' )


    end

end