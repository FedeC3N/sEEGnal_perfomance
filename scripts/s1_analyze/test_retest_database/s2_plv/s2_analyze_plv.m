%{

Evaluate differnces between sEEGnal and human experts.
 
For each subject, we estimate correlation, NRMSE, t-test sensor by sensor
and save the results.


@author: Fede

%}

clear
clc
restoredefaultpath

% Add paths
addpath('../shared/')

% Paths
config.path.clean_data = '../../../../databases/test_retest_database/derivatives';
config.path.results = '../../../../results/test_retest_database/plv';
if ~exist(config.path.results), mkdir(config.path.results),end

% Desired tasks
config.tasks = {'EO', 'EC'};

% To define later the pow matrix
config.complete_channel_labels = {'Fp1', 'Fpz', 'Fp2', 'F7', 'F3', 'Fz', 'F4', 'F8', 'FC5', 'FC1', 'FC2', 'FC6',...
    'T7', 'C3', 'Cz', 'C4', 'T8', 'CP5', 'CP1', 'CP2', 'CP6', 'P7', 'P3',...
    'Pz', 'P4', 'P8', 'POz', 'O1', 'O2', 'AF7', 'AF3', 'AF4', 'AF8', 'F5', 'F1', 'F2',...
    'F6', 'FC3', 'FCz', 'FC4', 'C5', 'C1', 'C2', 'C6', 'CP3', 'CP4', 'P5', 'P1', 'P2',...
    'P6', 'F9', 'PO3', 'PO4', 'F10', 'FT7', 'FT8', 'TP7', 'TP8', 'PO7', 'PO8', 'FT9',...
    'FT10', 'TPP9h', 'TPP10h', 'PO9', 'PO10', 'P9', 'P10', 'AFF1', 'AFz', 'AFF2', 'FFC5h',...
    'FFC3h', 'FFC4h', 'FFC6h', 'FCC5h', 'FCC3h', 'FCC4h', 'FCC6h', 'CCP5h', 'CCP3h', 'CCP4h',...
    'CCP6h', 'CPP5h', 'CPP3h', 'CPP4h', 'CPP6h', 'PPO1', 'PPO2', 'I1', 'Iz', 'I2', 'AFp3h', 'AFp4h',...
    'AFF5h', 'AFF6h', 'FFT7h', 'FFC1h', 'FFC2h', 'FFT8h', 'FTT9h', 'FTT7h', 'FCC1h', 'FCC2h', 'FTT8h',...
    'FTT10h', 'TTP7h', 'CCP1h', 'CCP2h', 'TTP8h', 'TPP7h', 'CPP1h', 'CPP2h', 'TPP8h', 'PPO9h', 'PPO5h',...
    'PPO6h', 'PPO10h', 'POO9h', 'POO3h', 'POO4h', 'POO10h', 'OI1h', 'OI2h'};

% Define the frequency bands
bands_info = struct('name',{'delta', 'theta','alpha','beta','gamma'},...
    'f_limits',{[ 2 4] , [4 8] , [8 12] , [12 20] , [20 45] });
config.bands_info = bands_info;

% Stats
% Output struct
results = [];
stats = []; % dimensions: [channels,recordings,bands,testers]
stats.NRMSE = nan(numel (config.complete_channel_labels),10,numel(bands_info),numel(config.tasks));
stats.rho = nan(numel (config.complete_channel_labels),10,numel(bands_info),numel(config.tasks));
stats.tstat = nan(numel (config.complete_channel_labels),10,numel(bands_info),numel(config.tasks));
stats.cohen_d = nan(numel (config.complete_channel_labels),10,numel(bands_info),numel(config.tasks));
stats.p = nan(numel (config.complete_channel_labels),10,numel(bands_info),numel(config.tasks));

for itask = 1 : numel(config.tasks)
    
    current_task = config.tasks{itask};

    % Differentiate EC and EO
    plv_dataset = read_plv_dataset(config, current_task);

    for iband = 1 : numel(bands_info)
        
        % For saving purposes
        current_band = bands_info(iband).name;
        
        for ichannel = 1 : numel(config.complete_channel_labels)
            
            % Get the current band matrix
            current_band_plv = nan(numel(config.complete_channel_labels),...
                size(plv_dataset,5),2);
            current_band_plv(:,:,1) = squeeze(plv_dataset(ichannel,:,iband,1,:));
            current_band_plv(:,:,2) = squeeze(plv_dataset(ichannel,:,iband,2,:));
            
            % Get the current PLV for a channel with the rest of the channels
            % Just interested in interconnectivity, not intra-connectivity
            current_band_plv(ichannel,:,:) = [];
            
            % For each subject, estimate the NMSE, corr, and t-test comparing
            % the two pre-processing pipelines
            for isubject = 1 : size(current_band_plv,2)
                
                % NMSE
                first = current_band_plv(:,isubject,1);
                second = current_band_plv(:,isubject,2);
                MSE = nanmean((first - second).^2);
                RMSE = sqrt(MSE);
                NRMSE = RMSE / 2; % To normalize we use the range, in PLV [-1 1]
                stats.NRMSE(ichannel,isubject,iband,itask) = NRMSE;
                
                % corr
                [rho,~] = corrcoef(first,second,'Rows','complete');
                stats.rho(ichannel,isubject,iband,itask) = rho(1,2);
                
                % t-test and  Effect Size (Cohen's d)
                diff = first - second;
                [~,p,~,t_stats] = ttest(diff);
                diff = abs(first - second);
                d = mean(diff)/std(diff);
                stats.tstat(ichannel,isubject,iband,itask) = t_stats.tstat;
                stats.cohen_d(ichannel,isubject,iband,itask) = d;
                stats.p(ichannel,isubject,iband,itask) = p;
                
            end
            
            
        end
       
    end
    
end

% Complete the information
results.stats = stats;
results.dimensions = 'channels x recordings x bands x task';
results.complete_channel_labels = config.complete_channel_labels;
results.bands_info = bands_info;
results.tasks = config.tasks;

% Save the file
outfile = sprintf('%s/plv_results.mat',config.path.results);
save(outfile,'-struct','results');


% Aux functions
function plv_dataset = read_plv_dataset(config,desired_task)

% Load the datset
dataset_path = sprintf('%s/sEEGnal/sEEGnal_dataset.mat',...
    config.path.clean_data);
dummy = load(dataset_path);

% Keep only the files of the condition
subjects = {dummy.dataset.sub};
subjects = unique(subjects);

% Create the output
% [channels x channels x bands_info x recordings x subjects]
plv_dataset = nan(numel(config.complete_channel_labels),numel(config.complete_channel_labels),...
    numel(config.bands_info),2,numel(subjects));

for isubject = 1 : numel(subjects)

    % Find the subject and 2 recordings of the current task
    dummy_subjects = {dummy.dataset.sub};
    dummy_tasks    = {dummy.dataset.task};
    dummy_tasks    = cellfun(@(x) x(2:end),dummy_tasks,'UniformOutput',false);
    subject_mask = ismember(dummy_subjects,subjects{isubject});
    task_mask = ismember(dummy_tasks,desired_task);
    dummy_mask = subject_mask & task_mask;
    dummy_index = find(dummy_mask);

    if numel(dummy_index) ~= 2
        continue
    end

    for irecording = 1 : numel(dummy_index)

        % Load pow
        current_dataset = dummy.dataset(dummy_index(irecording));
        plv_file = sprintf('../../../../%s/%s.mat',current_dataset.plv.path,...
            current_dataset.plv.file);
        plv = load(plv_file);

        for iband = 1 : numel(config.bands_info)

            % Add to the all matrix
            plv_dataset(:,:,iband,irecording,isubject) = plv.(config.bands_info(iband).name).plv;

        end


    end



end

end

