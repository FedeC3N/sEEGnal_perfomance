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

% HUMAN
fprintf('Humans\n')
data_average = nanmean(artifacts_human);
data_std = nanstd(artifacts_human);
fprintf('  Number: %.2f +- %.2f\n', ...
    data_average,data_std)

fprintf('sEEGnal\n')
data_average = nanmean(artifacts_sEEGnal);
data_std = nanstd(artifacts_sEEGnal);
fprintf('  Number: %.2f +- %.2f\n', ...
    data_average,data_std)


% Print statistical results
fprintf('\n\nSTATISTICAL RESULTS \n\n')


% Artifacts
[~,p,~,stats] = ttest2(artifacts_human,...
    artifacts_sEEGnal);
cohen_d = (2*stats.tstat)/sqrt(stats.df);

fprintf('Duration of Artifacts\n');
fprintf('   p-value: %.3f\n',p)
fprintf('   Cohen-d: %.3f\n',cohen_d)



% Read all the performance files
function artifacts_duration = read_performance_dataset_sEEGnal(config)

% Load the dataset
dataset_path = fullfile(config.path.derivatives,'sEEgnal','sEEGnal_dataset.mat');
dummy = load(dataset_path);

% Output structure
artifacts_duration = nan(numel(dummy.dataset),1);

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

    current_artifacts = [merged(:,1), merged(:,2) - merged(:,1)];
    current_total_duration = sum(current_artifacts,1);
    artifacts_duration(icurrent) = current_total_duration(2);    


end

end


function artifacts_duration = read_performance_dataset_human(config)

% Get the different testers (except sEEGnal)
testers = dir(sprintf('%s/*',config.path.derivatives));
testers = testers(3:end-1);
testers = {testers.name};

% Output structure
artifacts_duration = nan(20*numel(testers),1);
counter = 1;

for itester = 1 : numel(testers)

    % Load the dataset
    current_tester = testers{itester};
    dataset_path = fullfile(config.path.derivatives,current_tester,...
        sprintf('%s_dataset.mat',current_tester));
    dummy = load(dataset_path);

    for icurrent = 1 : numel(dummy.dataset)

        % Load performance
    current_dataset = dummy.dataset(icurrent);
    performance_file = fullfile(current_dataset.performance.path,...
        current_dataset.performance.file);
    current_performance = load(performance_file);

    % Transform to [begging, end]
    current_artifacts = [current_performance.artifacts.onset, ...
        current_performance.artifacts.onset + current_performance.artifacts.duration ];
    if numel(current_artifacts) == 0, continue, end
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

    current_artifacts = [merged(:,1), merged(:,2) - merged(:,1)];
    current_total_duration = sum(current_artifacts,1);
    artifacts_duration(counter)= current_total_duration(2);
    counter = counter + 1;
    end

end

end
