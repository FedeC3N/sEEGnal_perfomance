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
measures = {'times_seconds','memory_bytes','badchannels','artifacts','ICs_rejected'};
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
% Time
[~,p,~,stats] = ttest2(performance_human.times_seconds.badchannel_detection,...
    performance_sEEGnal.times_seconds.badchannel_detection);
cohen_d = (2*stats.tstat)/sqrt(stats.df);
fprintf('Time of Badchannels\n');
fprintf('   p-value: %.3f\n',p)
fprintf('   Cohen-d: %.3f\n',cohen_d)

[~,p,~,stats] = ttest2(performance_human.times_seconds.artifact_detection,...
    performance_sEEGnal.times_seconds.artifact_detection);
cohen_d = (2*stats.tstat)/sqrt(stats.df);
fprintf('Time of Artifact\n');
fprintf('   p-value: %.3f\n',p)
fprintf('   Cohen-d: %.3f\n',cohen_d)

% Memory
[~,p,~,stats] = ttest2(performance_human.memory_bytes.badchannel_detection,...
    performance_sEEGnal.memory_bytes.badchannel_detection);
cohen_d = (2*stats.tstat)/sqrt(stats.df);
fprintf('Memory of Badchannels\n');
fprintf('   p-value: %.3f\n',p)
fprintf('   Cohen-d: %.3f\n',cohen_d)

[~,p,~,stats] = ttest2(performance_human.memory_bytes.artifact_detection,...
    performance_sEEGnal.memory_bytes.artifact_detection);
cohen_d = (2*stats.tstat)/sqrt(stats.df);
fprintf('Memory of Artifact\n');
fprintf('   p-value: %.3f\n',p)
fprintf('   Cohen-d: %.3f\n',cohen_d)

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

fprintf('Number of ICs corrected\n');
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
performance.badchannels_description = cell(numel(dummy.dataset),124);

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
    
    % Store the badchannel description
    descriptions = current_performance.badchannels.description;
    performance.badchannels_description(icurrent,1:numel(descriptions)) = descriptions;
    
    % Count the number of rejected ICs and store it
    current_ICs_rejected = current_performance.ICs_rejected.IC;
    current_ICs_rejected = current_ICs_rejected(3:4);
    ICs_rejected_index = strcmp(current_ICs_rejected,'n/a');
    current_ICs_rejected = current_ICs_rejected(~ICs_rejected_index);
    dummy_count = strjoin(current_ICs_rejected,',');
    dummy_count = strsplit(dummy_count,',');
    performance.ICs_rejected(icurrent) = numel(dummy_count);
    
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
        for imodule = 1 : numel(modules)
            performance.times_seconds.(modules{imodule})(icurrent,itester) = current_performance.times_seconds.(modules{imodule});
            performance.memory_bytes.(modules{imodule})(icurrent,itester) = current_performance.memory_bytes.(modules{imodule}) * 0.000001;
        end
        
        % Store number of badchannels and artifacts
        performance.badchannels(icurrent,itester) = numel(current_performance.badchannels.label);
        performance.artifacts(icurrent,itester) = numel(current_performance.artifacts.label);
        
        % Count the number of rejected ICs and store it
        performance.ICs_rejected(icurrent,itester) = numel(current_performance.ICs_rejected.IC);
        
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

end
