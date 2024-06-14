clear
clc

dummy = dir("../../../aimind.etl.data/ETL_performance_AIMIND_dataset/raw/*cnt");
dummy = {dummy.name};

subjects = cellfun(@(x) x(1:9),dummy,'UniformOutput',false);
subjects = unique(subjects);


for isubject = 1:numel(subjects)
    
    current_subject = subjects{isubject};
    
    % Create the folder (and subfolder)
    new_folder = sprintf('./raw/prospective/eeg/%s/%s',...
        current_subject,current_subject);
    if ~exist(new_folder)
        mkdir(new_folder)
    end
    
    % List the files
    files_to_move = dir(sprintf('./raw/%s*',current_subject));
    files_to_move = {files_to_move.name};   
    
    % Move each file
    for ifile = 1:numel(files_to_move)
       
        current_file = files_to_move{ifile};
        
        old_file = sprintf('./raw/%s',current_file);
        new_file = sprintf('./raw/prospective/eeg/%s/%s/%s',...
            current_subject,current_subject,current_file);
        
        movefile(old_file,new_file);
        
        
    end
    
end
    