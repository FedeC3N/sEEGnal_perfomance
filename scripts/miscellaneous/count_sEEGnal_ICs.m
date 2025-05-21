%{
script to obtain clean recordings from the metadata

Created on Thu 16/04/2025

@author: Fede
%}
clear
clc
restoredefaultpath

% Paths
config.path.SOBI_files = './SOBI_files';

% Get the SOBI files
SOBI_files = dir('./SOBI_files/*.tsv');
SOBI_files = {SOBI_files.name};

% Matrix to save the count
n_components = zeros(1,numel(SOBI_files));

for ifile = 1 : numel(SOBI_files)

    % Read the current file
    current_file = sprintf('%s/%s',config.path.SOBI_files, SOBI_files{ifile});
    t = readtable(current_file, "FileType","text",'Delimiter', '\t');

    % Get the EOG and EKG rows
    index = ismember(table2cell(t(:,3)),{'eog','ecg'});
    t = t(index,4);

    % Convert to table and count the components
    t = table2cell(t);
    for i_comp = 1 : size(t,1)

        % Divide by comma
        current_comp = t{i_comp};

        % If no component, skip
        if ~strcmp(current_comp,'n/a')
            % Else, count
            current_comp = strsplit(current_comp,',');
            n_components(ifile) = n_components(ifile) + numel(current_comp);
        end


    end


end


fprintf('\n%.2f +- %.2f\n\n',mean(n_components),std(n_components))