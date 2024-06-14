% Each group's format is different.
% Unify group by group:
% ${country}_processed_${file_name}
clear
clc

Fs = 2000;
load('../channels_label.mat');


addpath('../../../../../SharedFunctions/functions/');

% Paths
path_outputs = '../../../data/metadata/original_format';
path_out = '../../../data/metadata/unified_format';

% File type
file_type = {'badchannels', 'annotations', 'ICs'};

% Finland
% badchannels and annotations
for itype = 1:2
    
    files = dir(sprintf('%s/Finland/%s/*_%s.tsv',path_outputs, file_type{itype},file_type{itype}));
    files = {files.name};
    
    for ifile = 1:numel(files)
        
        current_file = sprintf('%s/Finland/%s/%s',path_outputs, file_type{itype},files{ifile});
        
        % Create the new structure
        new_name = files{ifile};
        new_name = new_name(1:14);
        outfile = sprintf('%s/Finland_processed_%s',path_out,new_name);
        
        if exist([outfile '.mat'])
            new_struct = load(outfile);
        else
            new_struct = [];
        end
        
        % Read the data
        %[data, header, raw] = tsvread( file ) reads in text file with tab-seperated variables. default value for data is nan.
        %alternative input/output option is suppluying header strings
        %[col1, col2, col3, ..., header, raw] = tsvread( file, header1, header2, header3, ... )
        [~, ~, raw] = tsvread(current_file);
        
        switch itype
            case 1
                
                new_struct.(file_type{itype}).channels_label = channels_label;
                if size(raw,2) == 1
                    new_struct.(file_type{itype}).badchannels_labels = [];
                    new_struct.(file_type{itype}).badchannels_index = [];
                else
                    
                    current_badchannels_labels = raw(2:end,2);
                    current_badchannels_index = find(ismember(channels_label,current_badchannels_labels));
                    new_struct.(file_type{itype}).badchannels_labels = current_badchannels_labels;
                    new_struct.(file_type{itype}).badchannels_index = current_badchannels_index;
                end
                
            case 2
                
                current_artifact_type = raw(3:end,3);
                current_start_sample = raw(3:end,1);
                current_start_sample = cellfun(@str2double,current_start_sample);
                current_start_sample = int32(current_start_sample * Fs);
                duration = raw(3:end,2);
                duration = cellfun(@str2double,duration);
                current_end_sample = current_start_sample + int32(duration*Fs);
                
                new_struct.(file_type{itype}).artifact_type = current_artifact_type;
                new_struct.(file_type{itype}).start_sample = current_start_sample;
                new_struct.(file_type{itype}).end_sample = current_end_sample;
                
        end
        
        
        save(outfile,'-struct','new_struct')
        
        
    end
    
end

% ICs
files = dir(sprintf('%s/Finland/%s/*_%s.mat',path_outputs, file_type{3},file_type{3}));
files = {files.name};

for ifile = 1:numel(files)
    
    current_file = sprintf('%s/Finland/%s/%s',path_outputs, file_type{3},files{ifile});
    
    % Create the new structure
    new_name = files{ifile};
    new_name = new_name(1:14);
    outfile = sprintf('%s/Finland_processed_%s',path_out,new_name);
    
    if exist([outfile '.mat'])
        new_struct = load(outfile);
    else
        new_struct = [];
    end
    
    current_ICs = load(current_file);
    
    % Python's indexes start in 0, Matlab in 1. I have to add 1 to the
    % index.
    IC_index_matlab = current_ICs.IC_selection + 1;
    
    IC_excluded_mask = zeros(size(current_ICs.mixing_matrix,1),1);
    IC_excluded_mask(IC_index_matlab) = 1;
    IC_excluded_mask = boolean(IC_excluded_mask);
    IC_included_mask = ~IC_excluded_mask;
    
    new_struct.(file_type{3}).IC_included_mask = IC_included_mask;
    new_struct.(file_type{3}).mixing_matrix = current_ICs.mixing_matrix;
    new_struct.(file_type{3}).unmixing_matrix = current_ICs.unmixing_matrix;
    
    save(outfile,'-struct','new_struct')
    
end


