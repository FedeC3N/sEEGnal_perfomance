%{

Assess differences in artifact duration per recording

@author: Fede

%}

clear
clc
restoredefaultpath

% Add functions to the path
addpath('../shared/')

% Paths
config.path.derivatives = fullfile('..','..','..','..','databases',...
    'AI_Mind_database','derivatives');

% Read the performance
badchannels_label_human = read_performance_dataset_human(config);
badchannels_label_sEEGnal = read_performance_dataset_sEEGnal(config);

recall_all = nan(numel(badchannels_label_sEEGnal),1);
precision_all = nan(numel(badchannels_label_sEEGnal),1);
jaccard_all = nan(numel(badchannels_label_sEEGnal),1);

for ifile = 1 : numel(badchannels_label_sEEGnal)

    [recall,precision,jaccard] = overlap_metrics(badchannels_label_sEEGnal{ifile},...
        badchannels_label_human{ifile});

    recall_all(ifile) = recall;
    precision_all(ifile) = precision;
    jaccard_all(ifile) = jaccard;

end

recall_all = nanmean(recall_all);
precision_all = nanmean(precision_all);
jaccard_all = nanmean(jaccard_all);

% Mostrar resultados
fprintf('BADCHANNEL OVERLAP\n\n')
fprintf('Recall (%% benchmark cubierto): %.2f%%\n', 100*recall_all);
fprintf('Precision (%% detección válida): %.2f%%\n', 100*precision_all);
fprintf('Jaccard (acuerdo global): %.2f%%\n', 100*jaccard_all);


% Read all the performance files
function badchannels_labels = read_performance_dataset_sEEGnal(config)

% Load the dataset
dataset_path = fullfile(config.path.derivatives,'sEEgnal','sEEGnal_dataset.mat');
dummy = load(dataset_path);

% Output structure
badchannels_labels = cell(numel(dummy.dataset),1);

for icurrent = 1 : numel(dummy.dataset)

    % Load performance
    current_dataset = dummy.dataset(icurrent);
    performance_file = fullfile(current_dataset.performance.path,...
        current_dataset.performance.file);
    current_performance = load(performance_file);

    % Save the badhcannels 
    badchannels_labels(icurrent) = {current_performance.badchannels.label(:)};

end

end


function badchannels_labels = read_performance_dataset_human(config)

% Get the different testers (except sEEGnal)
testers = dir(sprintf('%s/*',config.path.derivatives));
testers = testers(3:end-1);
testers = {testers.name};

% Output structure
badchannels_labels = cell(20,1);

for icurrent = 1 : 20

    for itester = 1 : numel(testers)

        % Load the dataset
        current_tester = testers{itester};
        dataset_path = fullfile(config.path.derivatives,current_tester,...
            sprintf('%s_dataset.mat',current_tester));
        dummy = load(dataset_path);

        % Load performance
        current_dataset = dummy.dataset(icurrent);
        performance_file = fullfile(current_dataset.performance.path,...
            current_dataset.performance.file);
        current_performance = load(performance_file);

        % Save the badhcannels
        if itester == 1
            badchannels_labels(icurrent) = { current_performance.badchannels.label(:)};
        else
            badchannels_labels{icurrent} = ...
    [badchannels_labels{icurrent}; current_performance.badchannels.label(:)];
        end
    
    end

end

end


function [recall,precision,jaccard] = overlap_metrics(A,B)

% A = cell array a probar
% B = cell array benchmark

A_unique = unique(A);   % elimina duplicados de A
B_unique = unique(B);   % elimina duplicados de B

common = intersect(A,B);

n_common = numel(common);
n_A = numel(A);
n_B = numel(B);
n_union = numel(union(A,B));

% Recall
if n_B == 0
    recall = 0;
else
    recall = n_common / n_B;
end

% Precision
if n_A == 0
    precision = 0;
else
    precision = n_common / n_A;
end

% Jaccard
if n_union == 0
    jaccard = 1;   % ambos vacíos → acuerdo perfecto
else
    jaccard = n_common / n_union;
end
end
