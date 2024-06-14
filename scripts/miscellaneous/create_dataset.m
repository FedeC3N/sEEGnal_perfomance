clear
clc

% Paths
path = [];
path.metadata = '../../data/metadata/unified_format';
path.dataset = '../../data/metadata/dataset';
path.eeg = '../../data/raw/prospective/eeg';




files = dir(sprintf('%s/*.mat',path.metadata));
files = {files.name};

my_dataset = struct('eeg_path',[],'eeg_file',[],'metadata_path',[],...
    'metadata_file',[],'processed_by',[],'subject',[],'visit',[],...
    'letter',[],'task',[]);

for ifile = 1 :numel(files)
    
    % Parse the file
    current_file = files{ifile};
    pattern = '(\w*_processed)_([1-4]{1}-[0-9]{3})-([0-9]{1})-([A-Z]{1})_([1-4]{1}-[A-Z]{2}).mat';
    tokens = regexp(current_file,pattern,'tokens');
    
    % Create the information
    processed_by = tokens{1}{1};
    subject = tokens{1}{2};
    visit = tokens{1}{3};
    letter = tokens{1}{4};
    task = tokens{1}{5};
    eeg_path = sprintf('%s/%s-%s-%s/%s-%s-%s',path.eeg,...
        subject,visit,letter,subject,visit,letter);
    eeg_file = dir(sprintf('%s/%s-%s-%s/%s-%s-%s/*%s*',path.eeg,...
        subject,visit,letter,subject,visit,letter,task));
    eeg_file = eeg_file.name;
    metadata_path = path.metadata;
    metadata_file = sprintf('%s_%s-%s-%s_%s.mat',processed_by,subject,...
        visit,letter,task);
    
    
    % Store the information
    my_dataset(ifile).eeg_path = eeg_path;
    my_dataset(ifile).eeg_file = eeg_file;
    my_dataset(ifile).metadata_path = metadata_path;
    my_dataset(ifile).metadata_file = metadata_file;
    my_dataset(ifile).processed_by = processed_by;
    my_dataset(ifile).subject = subject;
    my_dataset(ifile).visit = visit;
    my_dataset(ifile).letter = letter;
    my_dataset(ifile).task = task;
    
end

outfile = sprintf('%s/dataset_all.mat',path.dataset);
save(outfile,'my_dataset')

