clear
clc
restoredefaultpath

% Select user
addpath('./.private')
user = validate_user();

%%%%%%%%
% Generate scripts

% Paths
config.path.scripts = sprintf('./scripts_%s',user);

% Copy the scripts folders
copyfile('./.private/common_scripts/', config.path.scripts)

% Create the user file
userfile_name = sprintf('%s/user.txt',config.path.scripts);
fid = fopen(userfile_name,'w');
fprintf(fid,'%s\n',user);
fclose(fid);


%%%%%%%%
% Generate dataset

% Paths
config.path.raw  = '../../data/raw';
config.path.dataset = '../../meta/dataset';
config.path.patt = '*.eeg';

% Create the output folder
if ~exist(config.path.dataset)
    mkdir(config.path.dataset)
end

% Lists the files in the EEG folder.
files = dir (sprintf('%s/%s',config.path.raw,config.path.patt) );
dataset    = struct ( 'file', [], 'subject', [], 'session',[],'task', []);
for ifile = 1 : length(files)
    
    % Extracts the information.
    file = files ( ifile ).name;
    parts = strsplit(file, '_');
    dummy = parts{1};
    dummy = strsplit(dummy,'-');
    subject = dummy{2};
    dummy = parts{2};
    dummy = strsplit(dummy,'-');
    session = dummy{2};
    dummy = parts{3};
    dummy = strsplit(dummy,'-');
    task = dummy{2};
    
    dataset(ifile).file = file;
    dataset(ifile).subject = subject;
    dataset(ifile).session = session;
    dataset(ifile).task = task;
    
end

% Saves the configuration structure.
outfile  = sprintf('%s/dataset_%s.mat',config.path.dataset,user);
save ( '-v6', outfile, 'dataset' )

%%%%%%%%
% Generate Fieldtrip metadata

% Adds the functions folders to the path.
addpath ( '../SharedFunctions/functions/');
addpath ( '../SharedFunctions/mne_silent/');
addpath ( '../SharedFunctions/functions_eep/');
addpath ( '../SharedFunctions/fieldtrip-20220626/');
ft_defaults

% Paths
config.path.raw  = '../../data/raw';
config.path.dataset = ['../../meta/dataset/dataset_' user '.mat'];
config.path.trl = ['../../meta/trl/' user];

% Action when the task has already been processed.
config.overwrite      = false;

% Defines the physiological channels (later will be relabeled).
config.physio.EOG     = { 'VEOGL' };
config.physio.EKG     = { 'CLAV' };
config.physio.EMG     = {};

% Sets the artifact detection parameters.
config.trialfun       = 'restingSegmentation';
config.segment        = 4;
config.padding        = 2;
config.addpadd        = true;
config.equal          = true;

% Creates the output folder, if needed.
if ~exist ( config.path.trl, 'dir' ), mkdir ( config.path.trl ); end

% Loads all the files.
dataset    = struct2array ( load ( config.path.dataset ) );

for ifile = 1 : numel(dataset)
    
    % Output information
    current_file = dataset(ifile).file;
    current_subject = dataset(ifile).subject;
    current_session = dataset(ifile).session;
    current_task = dataset(ifile).task;
    current_stage = [];
    fileinfo = [];
    artinfo = [];
    chaninfo = [];
    history = [];
    
    % Check if overwrite
    outfile   = sprintf ( '%s/sub-%s_ses-%s_task-%s_eeg.mat', config.path.trl, ...
        current_subject, current_session, current_task );
    fprintf('Current file %s\n', current_file);
    if exist(outfile) && ~config.overwrite
        fprintf('  Already exist. Do not overwrite\n\n');
        continue 
    end
    
    % fileinfo
    % 'file','subject','task','stage','begtime','endtime','header','event'
    fileinfo.file = current_file;
    fileinfo.subject = current_subject;
    fileinfo.session = current_session;
    fileinfo.task = current_task;
    fileinfo.stage = current_stage;
    
    % Read the header
    dummy_path_current_file = sprintf('%s/%s',config.path.raw,current_file);
    header = ft_read_header(dummy_path_current_file);
    fileinfo.header = header;
    
    % Look for events
    event = ft_read_event(dummy_path_current_file, 'header',header);
    fileinfo.event = event;
    
    % There is no trigger for begtime or endtime
    fileinfo.begtime = floor(1/header.Fs);
    fileinfo.endtime = floor(header.nSamples/header.Fs);
    
    % artifact
    % 'eog','muscle','jump','visual'
    artinfo.artifact.eog.artifact = zeros(0,2);
    artinfo.artifact.muscle.artifact = zeros(0,2);
    artinfo.artifact.jump.artifact = zeros(0,2);
    artinfo.artifact.visual.artifact = zeros(0,2);

    % chaninfo
    % 'bad'
    chaninfo.bad      = {''};
    
    % istory
    % 'step','date','history'
    history.last_step         = 'Creating Metadata';
    history.date         = datestr ( now );
    history.previous_steps      = {history.last_step};
    
    % Sets the task information.
    metadata          = [];
    metadata.subject  = current_subject;
    metadata.session  = current_session;
    metadata.task     = current_task;
    metadata.stage    = current_stage;
    metadata.fileinfo = fileinfo;
    metadata.chaninfo = chaninfo;
    metadata.history  = history;
    metadata.artinfo = artinfo;
    
    % Saves the output data.
    save ( '-v6', outfile, '-struct', 'metadata' );
    
end

