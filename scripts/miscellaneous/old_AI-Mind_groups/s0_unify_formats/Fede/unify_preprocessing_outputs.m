% Each group's format is different.
% Unify group by group:
% ${country}_processed_${file_name}
clear
clc

Fs = 2000;
load('../channels_labels.mat');


addpath('../../../../../SharedFunctions/functions/');

% Paths
path_outputs = '../../../data/metadata/original_format';
path_out = '../../../data/metadata/unified_format';

% File type
file_type = {'badchannels', 'annotations', 'ICs'};

% Fede
% badchannels, annotations and ICs
for itype = 1:2
    
    % There is a problem with read cell if the path is too long. I have to
    % move to the folder itself
    original_folder = pwd;
    cd(sprintf('%s/Fede/%s/',path_outputs, file_type{itype}));
    
    files = dir(sprintf('./*_%s.xlsx',file_type{itype}));
    files = {files.name};
    
    if isempty(files)
        files = dir(sprintf('./*_%s.xls',file_type{itype}));
        files = {files.name};
    end
    
    for ifile = 1:numel(files)
        
        current_file = files{ifile};
        
        % Create the new structure
        new_name = files{ifile};
        new_name = new_name(1:14);
        outfile = sprintf('../../../unified_format/Fede_processed_%s',new_name);
        
        if exist([outfile '.mat'])
            new_struct = load(outfile);
        else
            new_struct = [];
        end
        
        
        switch itype
            case 1
                
                % Read the data
                try
                    current_badchannels_labels = readcell(current_file,'Range','B1');
                    current_badchannels_labels = cellfun(@strsplit,current_badchannels_labels,'UniformOutput',false);
                    current_badchannels_labels = current_badchannels_labels{1};
                catch
                    current_badchannels_labels = [];
                end
                
                new_struct.(file_type{itype}).channel_labels = channel_labels;
                if isempty(current_badchannels_labels)
                    new_struct.(file_type{itype}).badchannels_labels = [];
                    new_struct.(file_type{itype}).badchannels_index = [];
                else
                    current_badchannels_index = find(ismember(channel_labels,current_badchannels_labels));
                    new_struct.(file_type{itype}).badchannels_labels = current_badchannels_labels';
                    new_struct.(file_type{itype}).badchannels_index = current_badchannels_index;
                end
                
                
            case 2
                
                raw = readcell(current_file);
                
                current_artifact_type = raw(2:end,1);
                current_start_sample = cell2mat(raw(2:end,2));
                current_end_sample = cell2mat(raw(2:end,3));
                
                new_struct.(file_type{itype}).artifact_type = current_artifact_type;
                new_struct.(file_type{itype}).start_sample = current_start_sample;
                new_struct.(file_type{itype}).end_sample = current_end_sample;
                
                
           
                
        end
        
        
        save(outfile,'-struct','new_struct')
        
    end
    
    cd(original_folder)
    
    
end

% ICs
files = dir(sprintf('%s/Fede/%s/*_%s.mat',path_outputs, file_type{3},file_type{3}));
files = {files.name};

for ifile = 1:numel(files)
    
    current_file = sprintf('%s/Fede/%s/%s',path_outputs, file_type{3},files{ifile});
    
    % Create the new structure
    new_name = files{ifile};
    new_name = new_name(1:14);    
    outfile = sprintf('%s/Fede_processed_%s',path_out,new_name);
    
    if exist([outfile '.mat'])
        new_struct = load(outfile);
    else
        new_struct = [];
    end
    
    current_ICs = load(current_file);
    
    IC_included_mask = current_ICs.IC_selection == 0;
    
    new_struct.(file_type{3}).IC_included_mask = IC_included_mask;
    new_struct.(file_type{3}).IC_types = current_ICs.IC_selection;
    new_struct.(file_type{3}).IC_types_description = current_ICs.IC_types;
    new_struct.(file_type{3}).mixing_matrix = current_ICs.mixing;
    new_struct.(file_type{3}).unmixing_matrix = current_ICs.unmixing;
    
    save(outfile,'-struct','new_struct')
    
end


