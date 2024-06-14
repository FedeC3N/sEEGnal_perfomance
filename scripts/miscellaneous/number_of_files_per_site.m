clear
clc

% Paths
path = [];
path.metadata = '../../data/metadata/unified_format';
path.dataset = '../../data/metadata/dataset';

% Files (based on dataset)
load(sprintf('%s/dataset_all.mat',path.dataset));

% sites
sites = {'Fede','Finland','IClabel','Italy'};

% number of files
files_per_site = nan(1,numel(sites));

% Go through sites
for isite = 1:numel(sites)
    
    % Keep only files preprocessed by the current site
    current_site = sites{isite};
    dummy_processed_str = sprintf('%s_processed',current_site);
    files_included_mask = strcmp({my_dataset.processed_by},dummy_processed_str);
   
    files_per_site(isite) = sum(files_included_mask);
    
end

% unique subjects per site
subjects_per_site = nan(1,numel(sites));

% Go through sites
for isite = 1:numel(sites)
    
    % Keep only files preprocessed by the current site
    current_site = sites{isite};
    dummy_processed_str = sprintf('%s_processed',current_site);
    files_included_mask = strcmp({my_dataset.processed_by},dummy_processed_str);
    
    subjects = {my_dataset(files_included_mask).subject};
    subjects = unique(subjects)';
    
   
    subjects_per_site(isite) = numel(subjects);
    
end

% Find the non-processed subjects 
% Processed
% Keep only files preprocessed by the current site
files_included_mask = strcmp({my_dataset.processed_by},'Fede_processed');
correct_subset = my_dataset(files_included_mask);
correct_subset_subjects = {correct_subset.subject}';
correct_subset_subjects = unique(correct_subset_subjects);

% Compare Finland
files_included_mask = strcmp({my_dataset.processed_by},'Finland_processed');
Finland_subset = my_dataset(files_included_mask);
Finland_subset_subjects = {Finland_subset.subject}';
Finland_subset_subjects = unique(Finland_subset_subjects);

missing_index = ~ismember(correct_subset_subjects,Finland_subset_subjects);
Finland_missing_subjects = correct_subset_subjects(missing_index);


% Compare IClabel
files_included_mask = strcmp({my_dataset.processed_by},'IClabel_processed');
IClabel_subset = my_dataset(files_included_mask);
IClabel_subset_subjects = {IClabel_subset.subject}';
IClabel_subset_subjects = unique(IClabel_subset_subjects);

missing_index = ~ismember(correct_subset_subjects,IClabel_subset_subjects);
IClabel_missing_subjects = correct_subset_subjects(missing_index);


% Find the non-processed recordings 
% Processed
% Keep only files preprocessed by the current site
files_included_mask = strcmp({my_dataset.processed_by},'Fede_processed');
correct_subset = my_dataset(files_included_mask);
correct_subset_files = {correct_subset.eeg_file}';

% Compare Finland
files_included_mask = strcmp({my_dataset.processed_by},'Finland_processed');
Finland_subset = my_dataset(files_included_mask);
Finland_subset_files = {Finland_subset.eeg_file}';

missing_index = ~ismember(correct_subset_files,Finland_subset_files);
Finland_missing_files = Finland_subset_files(missing_index);


% Compare IClabel
files_included_mask = strcmp({my_dataset.processed_by},'IClabel_processed');
IClabel_subset = my_dataset(files_included_mask);
IClabel_subset_files = {IClabel_subset.eeg_file}';

missing_index = ~ismember(correct_subset_files,IClabel_subset_files);
IClabel_missing_files = IClabel_subset_files(missing_index);
