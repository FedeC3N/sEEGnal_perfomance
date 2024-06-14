% Each group's format is different.
% Unify group by group:
% ${country}_processed_${file_name}
clear
clc

Fs = 2000;
load('../channels_label_Italy.mat');


addpath('../../../../../SharedFunctions/functions/');

% Paths
path_outputs = '../../../data/metadata/original_format';
path_out = '../../../data/metadata/unified_format';

% File type
file_type = {'badchannels', 'annotations', 'ICs'};

% badchannels, annotations and ICs
for itype = 1:3
    
    % There is a problem with read cell if the path is too long. I have to
    % move to the folder itself
    original_folder = pwd;
    cd(sprintf('%s/Italy/%s/',path_outputs, file_type{itype}));
    
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
        outfile = sprintf('../../../unified_format/Italy_processed_%s',new_name);
        
        if exist([outfile '.mat'])
            new_struct = load(outfile);
        else
            new_struct = [];
        end
        
        
        switch itype
            case 1
                
                % Read the data
                current_badchannels_labels = readcell(current_file,'Range','E2');
                
                new_struct.(file_type{itype}).channel_labels = channels_label;
                if strcmp(current_badchannels_labels,'n/a')
                    new_struct.(file_type{itype}).badchannels_labels = [];
                    new_struct.(file_type{itype}).badchannels_index = [];
                else
                    current_badchannels_index = find(ismember(channels_label,current_badchannels_labels));
                    new_struct.(file_type{itype}).badchannels_labels = current_badchannels_labels';
                    new_struct.(file_type{itype}).badchannels_index = current_badchannels_index;
                end
                
                
            case 2
                
                raw = readcell(current_file);
                
                current_artifact_type = raw(2:end,4);
                current_start_sample = raw(2:end,3);
                current_start_sample = cell2mat(current_start_sample);
                current_start_sample = current_start_sample + 1;
                duration = 2;
                current_end_sample = current_start_sample + duration*Fs;
                current_end_sample = current_end_sample + 1;
                
                new_struct.(file_type{itype}).artifact_type = current_artifact_type;
                new_struct.(file_type{itype}).start_sample = current_start_sample;
                new_struct.(file_type{itype}).end_sample = current_end_sample;
                
                
            case 3
                
                current_ICs = [];
                dummy_subject = current_file(1:14);
                
                
                current_ICs.IC_selection = readcell([dummy_subject '_ICs.xls'],...
                    'Range','E:E');
                current_ICs.IC_selection = cell2mat(current_ICs.IC_selection(2:end));
                current_ICs.mixing_matrix = readmatrix([dummy_subject '_ICAmixing.xls'],...
                    'Range', 'B2:DU125');
                current_ICs.unmixing_matrix = readmatrix([dummy_subject '_ICAunmixing.xls'],...
                    'Range', 'B2:DU125');
                
                IC_included_mask = ones(size(current_ICs.mixing_matrix,1),1);
                IC_included_mask(current_ICs.IC_selection) = 0;
                
                new_struct.(file_type{itype}).IC_selection = current_ICs.IC_selection;
                new_struct.(file_type{itype}).mixing_matrix = current_ICs.mixing_matrix;
                new_struct.(file_type{itype}).unmixing_matrix = current_ICs.unmixing_matrix;
                new_struct.(file_type{itype}).IC_included_mask = IC_included_mask;
                
                
        end
        
        
        save(outfile,'-struct','new_struct')
        
    end
    
    cd(original_folder)
    
    
end

