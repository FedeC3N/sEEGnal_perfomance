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

% Matrix to store the number of badchannels
badchannels_all = nan(50,numel(sites),2);

% Go through sites
for isite = 1:numel(sites)
    
    % Keep only files preprocessed by the current site
    current_site = sites{isite};
    dummy_processed_str = sprintf('%s_processed',current_site);
    files_site_mask = strcmp({my_dataset.processed_by},dummy_processed_str);
    
    % Keep only EO
    files_1EO_mask = strcmp({my_dataset.task},'1-EO');
    files_3EO_mask = strcmp({my_dataset.task},'3-EO');
    
    files_included_mask = files_site_mask & (files_1EO_mask | files_3EO_mask);
    current_subset = my_dataset(files_included_mask);
    
    % For each file
    for ifile = 1:numel(current_subset)
        
        % Load the current metadata
        current_metadata = current_subset(ifile);
        
        % Load the metadata and extract the information
        metadata = load(sprintf('%s/%s',current_metadata.metadata_path,...
            current_metadata.metadata_file));
        
        % Save the number of badchannels
        if isfield(metadata,'badchannels')
            badchannels_all(ifile,isite,1) = numel(metadata.badchannels.badchannels_index);
        else
            badchannels_all(ifile,isite,1) = 0;
        end
    end
    
    % Keep only EC
    files_2EC_mask = strcmp({my_dataset.task},'2-EC');
    files_4EC_mask = strcmp({my_dataset.task},'4-EC');
    
    files_included_mask = files_site_mask & (files_2EC_mask | files_4EC_mask);
    current_subset = my_dataset(files_included_mask);
    
    % For each file
    for ifile = 1:numel(current_subset)
        
        % Load the current metadata
        current_metadata = current_subset(ifile);
        
        % Load the metadata and extract the information
        metadata = load(sprintf('%s/%s',current_metadata.metadata_path,...
            current_metadata.metadata_file));
        
        % Save the number of badchannels
        if isfield(metadata,'badchannels')
            badchannels_all(ifile,isite,2) = numel(metadata.badchannels.badchannels_index);
        else
            badchannels_all(ifile,isite,2) = 0;
        end
    end
    
end

mean_badchannels = mean(badchannels_all,1,'omitnan');
std_badchannels = std(badchannels_all,1,'omitnan');

% Get the difference in number of channels among sites
badchannels_EO = badchannels_all(:,:,1);
dummy = reshape(badchannels_EO,size(badchannels_EO,1),1,size(badchannels_EO,2));
difference_EO = squeeze(mean(abs(badchannels_EO - dummy),1,'omitnan'));
difference_EO(difference_EO == 0) = nan;
mean_difference_EO = mean(difference_EO,1,'omitnan');

badchannels_EC = badchannels_all(:,:,2);
dummy = reshape(badchannels_EC,size(badchannels_EC,1),1,size(badchannels_EC,2));
difference_EC = squeeze(mean(abs(badchannels_EC - dummy),1,'omitnan'));
difference_EC(difference_EC == 0) = nan;
mean_difference_EC = mean(difference_EC,1,'omitnan');

