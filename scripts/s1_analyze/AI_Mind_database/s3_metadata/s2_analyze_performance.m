%{

Assess differences in: 
 
   {'times_seconds','memory_bytes','badchannels',
          'artifacts','ICs_rejected'};

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
performance_human = read_performance_dataset_human(config);
performance_sEEGnal = read_performance_dataset_sEEGnal(config);

% Analyze
modules = {'badchannel_detection', 'artifact_detection'  };
measures = {'times_seconds','memory_bytes','badchannels','artifacts','ICs_rejected','rejected_variance'};
for imeasure = 1 : numel(measures)
    
    % Current measure
    current_measure = measures{imeasure};
    fprintf('%s\n',current_measure);
    
    % Print
    % HUMAN
    fprintf('Humans\n')
    if imeasure < 3
        
        for imodule = 1 : numel(modules)
            current_module = modules{imodule};
            data_average = nanmean(performance_human.(current_measure).(current_module));
            data_std = nanstd(performance_human.(current_measure).(current_module));
            fprintf('  %s: %.2f +- %.2f\n', ...
                modules{imodule}, data_average,data_std)
        end
    else
        data_average = nanmean(performance_human.(current_measure));
        data_std = nanstd(performance_human.(current_measure));
        fprintf('  Number: %.2f +- %.2f\n', ...
            data_average,data_std)
    end
    fprintf('\n')
    
    % sEEGnal
    fprintf('sEEGnal\n')
    if imeasure < 3
        
        for imodule = 1 : numel(modules)
            current_module = modules{imodule};
            data_average = nanmean(performance_sEEGnal.(current_measure).(current_module));
            data_std = nanstd(performance_sEEGnal.(current_measure).(current_module));
            fprintf('  %s: %.2f +- %.2f\n', ...
                modules{imodule}, data_average,data_std)
        end
    else
        data_average = nanmean(performance_sEEGnal.(current_measure));
        data_std = nanstd(performance_sEEGnal.(current_measure));
        fprintf('  Number: %.2f +- %.2f\n', ...
            data_average,data_std)
    end
    fprintf('\n\n\n')
    
end

% Print statistical results
fprintf('\n\nSTATISTICAL RESULTS \n\n')

% T-test and effect size and print
% Badchannels
[~,p,~,stats] = ttest2(performance_human.badchannels,...
    performance_sEEGnal.badchannels);
cohen_d = (2*stats.tstat)/sqrt(stats.df);

fprintf('Number of Badchannels\n');
fprintf('   p-value: %.3f\n',p)
fprintf('   Cohen-d: %.3f\n',cohen_d)

% Artifacts
[~,p,~,stats] = ttest2(performance_human.artifacts,...
    performance_sEEGnal.artifacts);
cohen_d = (2*stats.tstat)/sqrt(stats.df);

fprintf('Number of Artifacts\n');
fprintf('   p-value: %.3f\n',p)
fprintf('   Cohen-d: %.3f\n',cohen_d)

% ICs
[~,p,~,stats] = ttest2(performance_human.ICs_rejected,...
    performance_sEEGnal.ICs_rejected);
cohen_d = (2*stats.tstat)/sqrt(stats.df);

fprintf('Number of ICs rejected\n');
fprintf('   p-value: %.3f\n',p)
fprintf('   Cohen-d: %.3f\n',cohen_d)

% ICs corrected
[~,p,~,stats] = ttest2(performance_human.ICs_rejected,...
    performance_sEEGnal.ICs_rejected_corrected);
cohen_d = (2*stats.tstat)/sqrt(stats.df);

fprintf('Number of ICs rejected corrected\n');
fprintf('   p-value: %.3f\n',p)
fprintf('   Cohen-d: %.3f\n',cohen_d)

% Variance rejected
[~,p,~,stats] = ttest2(performance_human.rejected_variance,...
    performance_sEEGnal.rejected_variance);
cohen_d = (2*stats.tstat)/sqrt(stats.df);

fprintf('Number of rejected variance\n');
fprintf('   p-value: %.3f\n',p)
fprintf('   Cohen-d: %.3f\n',cohen_d)


% Read all the performance files
function performance = read_performance_dataset_sEEGnal(config)

% Load the dataset
dataset_path = fullfile(config.path.derivatives,'sEEgnal','sEEGnal_dataset.mat');
dummy = load(dataset_path);

% Output structure
modules = {'badchannel_detection', 'artifact_detection'  };
performance = [];
for imodule = 1 : numel(modules)
    performance.times_seconds .(modules{imodule})= nan(numel(dummy.dataset),1);
    performance.memory_bytes.(modules{imodule}) = nan(numel(dummy.dataset),1);
end
performance.badchannels = nan(numel(dummy.dataset),1);
performance.artifacts = nan(numel(dummy.dataset),1);
performance.ICs_rejected = nan(numel(dummy.dataset),1);
performance.rejected_variance = nan(numel(dummy.dataset),1);

for icurrent = 1 : numel(dummy.dataset)
    
    % Load performance
    current_dataset = dummy.dataset(icurrent);
    performance_file = fullfile(current_dataset.performance.path,...
        current_dataset.performance.file);
    current_performance = load(performance_file);
    
    % Store times and bytes
    for imodule = 1 : numel(modules)
        performance.times_seconds.(modules{imodule})(icurrent) = current_performance.times_seconds.(modules{imodule});
        performance.memory_bytes.(modules{imodule})(icurrent) = current_performance.memory_bytes.(modules{imodule}) * 0.000001;
    end
    
    % Store number of badchannels and artifacts
    performance.badchannels(icurrent) = numel(current_performance.badchannels.label);
    performance.artifacts(icurrent) = numel(current_performance.artifacts.label);
    
    % Count the number of rejected ICs and store it
    current_ICs_rejected = current_performance.ICs_rejected.IC;
    current_ICs_rejected = current_ICs_rejected(2:end-1);
    ICs_rejected_index = strcmp(current_ICs_rejected,'n/a');
    current_ICs_rejected = current_ICs_rejected(~ICs_rejected_index);
    dummy_count = strjoin(current_ICs_rejected,',');
    dummy_count = strsplit(dummy_count,',');
    performance.ICs_rejected(icurrent) = numel(dummy_count);

    % Count the number of rejected ICs corrected to be compared to humans
    current_ICs_rejected_corrected = current_performance.ICs_rejected.IC;
    current_ICs_rejected_corrected = current_ICs_rejected_corrected(3:4);
    ICs_rejected_index = strcmp(current_ICs_rejected_corrected,'n/a');
    current_ICs_rejected_corrected = current_ICs_rejected_corrected(~ICs_rejected_index);
    dummy_count = strjoin(current_ICs_rejected_corrected,',');
    dummy_count = strsplit(dummy_count,',');
    performance.ICs_rejected_corrected(icurrent) = numel(dummy_count);

    % Sum the variance rejected and store it
    current_ICs_rejected = current_performance.ICs_rejected.IC;
    current_ICs_rejected = current_ICs_rejected(2:end-1);
    ICs_rejected_index = strcmp(current_ICs_rejected,'n/a');
    current_ICs_rejected = current_ICs_rejected(~ICs_rejected_index);
    current_ICs_rejected = strjoin(current_ICs_rejected,',');
    current_ICs_rejected = strsplit(current_ICs_rejected,',');
    rejected_index = cellfun(@(x) x(3:end),current_ICs_rejected,'UniformOutput',false);
    rejected_index = cellfun(@str2double,rejected_index,'UniformOutput',false);
    rejected_index = cell2mat(rejected_index);
    rejected_variance = sum(current_performance.ICs_rejected.explained_variance(rejected_index));
    performance.rejected_variance(icurrent) = rejected_variance;

end

end


function performance = read_performance_dataset_human(config)

% Get the different testers (except sEEGnal)
testers = dir(sprintf('%s/*',config.path.derivatives));
testers = testers(3:end-1);
testers = {testers.name};

% Output structure
modules = {'badchannel_detection', 'artifact_detection'  };
performance = [];
for imodule = 1 : numel(modules)
    performance.times_seconds .(modules{imodule})= nan(20,numel(testers));
    performance.memory_bytes.(modules{imodule}) = nan(20,numel(testers));
end
performance.badchannels = nan(20,numel(testers));
performance.artifacts = nan(20,numel(testers));
performance.ICs_rejected = nan(20,numel(testers));
performance.rejected_variance = nan(20,numel(testers));

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
        
        % Store times and bytes
        % for imodule = 1 : numel(modules)
            % performance.times_seconds.(modules{imodule})(icurrent,itester) = current_performance.times_seconds.(modules{imodule});
            % performance.memory_bytes.(modules{imodule})(icurrent,itester) = current_performance.memory_bytes.(modules{imodule}) * 0.000001;
        % end
        
        % Store number of badchannels and artifacts
        performance.badchannels(icurrent,itester) = numel(current_performance.badchannels.label);
        performance.artifacts(icurrent,itester) = numel(current_performance.artifacts.label);
        
        % Count the number of rejected ICs and store it
        performance.ICs_rejected(icurrent,itester) = numel(current_performance.ICs_rejected.IC);
        
        % Sum the variance rejected and store it
        rejected_index = current_performance.ICs_rejected.IC;
        rejected_variance = sum(current_performance.ICs_rejected.explained_variance(rejected_index));
        performance.rejected_variance(icurrent,itester) = rejected_variance;

    end
    
end

% Average
for imodule = 1 : numel(modules)
    performance.times_seconds.(modules{imodule}) = nanmean(performance.times_seconds.(modules{imodule}),2);
    performance.memory_bytes.(modules{imodule}) = nanmean(performance.memory_bytes.(modules{imodule}),2);
end
performance.badchannels = nanmean(performance.badchannels,2);
performance.artifacts = nanmean(performance.artifacts,2);
performance.ICs_rejected = nanmean(performance.ICs_rejected,2);
performance.rejected_variance = nanmean(performance.rejected_variance,2);

end
