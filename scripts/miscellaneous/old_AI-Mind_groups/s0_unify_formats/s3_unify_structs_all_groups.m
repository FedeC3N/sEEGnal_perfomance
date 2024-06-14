% Each group's format is different.
% Unify group by group:
% ${country}_processed_${file_name}
clear
clc

Fs = 2000;
load('./channels_labels.mat');


addpath('../../../../SharedFunctions/functions/');

% Paths
path_outputs = '../../data/original_format';
path_out = '../../data/unified_format';

% File type
file_type = {'badchannels', 'annotations', 'ICs'};

%% Finland
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
                
                new_struct.(file_type{itype}).channel_labels = channel_labels;
                if size(raw,2) == 1
                    new_struct.(file_type{itype}).badchannels_labels = [];
                    new_struct.(file_type{itype}).badchannels_index = [];
                else
                    
                    current_badchannels_labels = raw(2:end,2);
                    current_badchannels_index = find(ismember(channel_labels,current_badchannels_labels));
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
    
    new_struct.(file_type{3}).IC_selection = current_ICs.IC_selection';
    new_struct.(file_type{3}).mixing_matrix = current_ICs.mixing_matrix;
    new_struct.(file_type{3}).unmixing_matrix = current_ICs.unmixing_matrix;
    
    save(outfile,'-struct','new_struct')
    
end


%% Italy
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
        outfile = sprintf('../../../Italy_processed_%s',new_name);
        
        if exist([outfile '.mat'])
            new_struct = load(outfile);
        else
            new_struct = [];
        end
        
        
        switch itype
            case 1
                
                % Read the data
                current_badchannels_labels = readcell(current_file,'Range','E2');
                
                new_struct.(file_type{itype}).channel_labels = channel_labels;
                if strcmp(current_badchannels_labels,'n/a')
                    new_struct.(file_type{itype}).badchannels_labels = [];
                    new_struct.(file_type{itype}).badchannels_index = [];
                else
                    current_badchannels_index = find(ismember(channel_labels,current_badchannels_labels));
                    new_struct.(file_type{itype}).badchannels_labels = current_badchannels_labels';
                    new_struct.(file_type{itype}).badchannels_index = current_badchannels_index;
                end
                
                
            case 2
                
                raw = readcell(current_file);
                
                current_artifact_type = raw(2:end,4);
                current_start_sample = raw(2:end,3);
                current_start_sample = cell2mat(current_start_sample);
                duration = 2;
                current_end_sample = current_start_sample + duration*Fs;
                
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
                
                new_struct.(file_type{itype}).IC_selection = current_ICs.IC_selection;
                new_struct.(file_type{itype}).mixing_matrix = current_ICs.mixing_matrix;
                new_struct.(file_type{itype}).unmixing_matrix = current_ICs.unmixing_matrix;
                
                
                
        end
        
        
        save(outfile,'-struct','new_struct')
        
    end
    
    cd(original_folder)
    
    
end

%% Fede
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
    
    new_struct.(file_type{3}).IC_selection = current_ICs.IC_selection;
    new_struct.(file_type{3}).IC_types = current_ICs.IC_types;
    new_struct.(file_type{3}).mixing_matrix = current_ICs.mixing;
    new_struct.(file_type{3}).unmixing_matrix = current_ICs.unmixing;
    
    save(outfile,'-struct','new_struct')
    
end


%% Ricardo
% badchannels, annotations and ICs
for itype = 1:2
    
    % There is a problem with read cell if the path is too long. I have to
    % move to the folder itself
    original_folder = pwd;
    cd(sprintf('%s/Ricardo/%s/',path_outputs, file_type{itype}));
    
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
        outfile = sprintf('../../../unified_format/Ricardo_processed_%s',new_name);
        
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



