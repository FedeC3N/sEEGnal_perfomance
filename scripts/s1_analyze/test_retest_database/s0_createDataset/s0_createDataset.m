%{

Create a dataset with the files to process

@author: Fede

%}

clear
clc
restoredefaultpath

% Paths
config.path.clean_data = '../../../../databases/test_retest_database/derivatives';

% Get the different testers
testers = dir(sprintf('%s/*',config.path.clean_data));
testers = testers(3:end);

for itester = 1 : numel(testers)
    
    % Structs for the files
    dataset = struct('origin',[],'sub',[],'ses',[],'task',[],'path',[],'file',[]);
    dataset_index = 1;
    
    % Get the files processed
    current_tester = testers(itester);
    subjects = dir(sprintf('%s/%s/clean/sub*',current_tester.folder, current_tester.name));
    
    for isubject = 1 : numel(subjects)

        % Find the session in the folder
        current_subject = subjects(isubject);
        current_session = dir(sprintf('%s/%s/ses*',current_subject.folder, current_subject.name));

        % Find the files to work with
        files = dir(sprintf('%s/%s/%s/eeg/*set',current_subject.folder, current_subject.name, current_session.name));
        
        for ifile = 1 : numel(files)
            
            % Get the number of session for the current subject
            current_file = files(ifile);
            
            % Get the metadata
            expression = 'sub-(\w*)_ses-(\w*)_task-(\w*)_desc-(\w*)_clean_eeg.set';
            dummy = regexp(current_file.name,expression,'tokens');
            
            % For the EEG experts, use 'eeg_expert'
            dataset(dataset_index).origin = dummy{1}{4};
            dataset(dataset_index).sub = dummy{1}{1};
            dataset(dataset_index).ses = dummy{1}{2};
            dataset(dataset_index).task = dummy{1}{3};
            dataset(dataset_index).file = current_file.name;
            dataset(dataset_index).path = fullfile('databases',...
                'test_retest_database','derivatives',current_tester.name,'clean',...
                sprintf('sub-%s',dataset(dataset_index).sub),...
                sprintf('ses-%s',dataset(dataset_index).ses),...
                'eeg');
            
            % Update the struct index
            dataset_index = dataset_index + 1;
            
            
        end
        
        
    end
    
    % Save the struct
    outfile = sprintf('%s/%s/%s_dataset.mat',current_tester.folder,...
        current_tester.name,current_tester.name);
    save(outfile,'dataset')
    
end






