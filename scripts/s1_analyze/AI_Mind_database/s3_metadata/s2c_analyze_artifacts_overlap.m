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
artifacts_human = read_performance_dataset_human(config);
artifacts_sEEGnal = read_performance_dataset_sEEGnal(config);

recall_all = nan(numel(artifacts_sEEGnal),1);
precision_all = nan(numel(artifacts_sEEGnal),1);
jaccard_all = nan(numel(artifacts_sEEGnal),1);


for ifile = 1 : numel(artifacts_sEEGnal)

    [recall,precision,jaccard] = overlap_metrics(artifacts_sEEGnal(ifile).definition,...
        artifacts_human(ifile).definition);

    recall_all(ifile) = recall;
    precision_all(ifile) = precision;
    jaccard_all(ifile) = jaccard;

end


recall_all = nanmean(recall_all);
precision_all = nanmean(precision_all);
jaccard_all = nanmean(jaccard_all);

% Mostrar resultados
fprintf('ARTIFACT OVERLAP\n\n')
fprintf('Recall (%% benchmark cubierto): %.2f%%\n', 100*recall_all);
fprintf('Precision (%% detección válida): %.2f%%\n', 100*precision_all);
fprintf('Jaccard (acuerdo global): %.2f%%\n', 100*jaccard_all);


% Read all the performance files
function artifacts = read_performance_dataset_sEEGnal(config)

% Load the dataset
dataset_path = fullfile(config.path.derivatives,'sEEgnal','sEEGnal_dataset.mat');
dummy = load(dataset_path);

% Output structure
artifacts = [];

for icurrent = 1 : numel(dummy.dataset)

    % Load performance
    current_dataset = dummy.dataset(icurrent);
    performance_file = fullfile(current_dataset.performance.path,...
        current_dataset.performance.file);
    current_performance = load(performance_file);

    % Transform to [begging, end]
    current_artifacts = [current_performance.artifacts.onset, ...
        current_performance.artifacts.onset + current_performance.artifacts.duration ];
    current_artifacts = sortrows(current_artifacts,1);

    % Merged overlapping
    merged = current_artifacts(1,:);

    for i = 2:size(current_artifacts,1)

        current = current_artifacts(i,:);
        last = merged(end,:);

        % Si hay solapamiento (incluye eventos que se tocan)
        if current(1) <= last(2)
            merged(end,2) = max(last(2), current(2));
        else
            merged = [merged; current];
        end
    end

    % ==============================
    % VOLVER A [onset, duration]
    % ==============================

    artifacts(icurrent).definition = merged;


end

end


function artifacts = read_performance_dataset_human(config)

% Get the different testers (except sEEGnal)
testers = dir(sprintf('%s/*',config.path.derivatives));
testers = testers(3:end-1);
testers = {testers.name};

% Output structure
artifacts = [];

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

        % Transform to [begging, end]
        if itester == 1
            current_artifacts = [current_performance.artifacts.onset, ...
                current_performance.artifacts.onset + current_performance.artifacts.duration ];
            if numel(current_artifacts) == 0

                current_artifacts = [-5,0.5];

            end
        else
            dummy = [current_performance.artifacts.onset, ...
                current_performance.artifacts.onset + current_performance.artifacts.duration ];
            if numel(dummy) == 0

                continue

            end
            current_artifacts = cat(1,current_artifacts,dummy);
        end


    end

    % Merged overlapping
    current_artifacts = sortrows(current_artifacts,1);
    merged = current_artifacts(1,:);

    for i = 2:size(current_artifacts,1)

        current = current_artifacts(i,:);
        last = merged(end,:);

        % Si hay solapamiento (incluye eventos que se tocan)
        if current(1) <= last(2)
            merged(end,2) = max(last(2), current(2));
        else
            merged = [merged; current];
        end
    end

    artifacts(icurrent).definition = merged;

end

end

function [recall,precision,jaccard] = overlap_metrics(A,B)

% ==========================
% INPUT
% ==========================
% A = detectado  [inicio, fin]
% B = benchmark  [inicio, fin]

A = sortrows(A,1);
B = sortrows(B,1);

% ==========================
% INTERSECCIÓN
% ==========================

i = 1;
j = 1;
overlaps = [];

while i <= size(A,1) && j <= size(B,1)

    start_overlap = max(A(i,1), B(j,1));
    end_overlap   = min(A(i,2), B(j,2));

    if start_overlap < end_overlap
        overlaps = [overlaps; start_overlap, end_overlap];
    end

    if A(i,2) < B(j,2)
        i = i + 1;
    else
        j = j + 1;
    end
end

% ==========================
% MÉTRICAS
% ==========================
if isempty(overlaps)
    intersection = 0;
else
    intersection = sum(overlaps(:,2) - overlaps(:,1));
end

total_A = sum(A(:,2) - A(:,1));
total_B = sum(B(:,2) - B(:,1));

% Evitar divisiones por cero
if total_B == 0
    recall = 0;
else
    recall = intersection / total_B;
end

if total_A == 0
    precision = 0;
else
    precision = intersection / total_A;
end

union = total_A + total_B - intersection;

if union == 0
    jaccard = 0;
else
    jaccard = intersection / union;
end

end
