clear
clc
restoredefaultpath

% Paths
config.path.clean_data = '../../databases/LEMON_database/derivatives/lemon/clean';

% Get the different folders
folders = dir(sprintf('%s/sub*',config.path.clean_data));

for ifolder = 1 : numel(folders)
   
    files = dir(sprintf('%s/%s/sub*',folders(ifolder).folder,...
        folders(ifolder).name));
    
    
    for ifile = 1 : numel(files)
        
        destination_folder = sprintf('%s/ses-1/eeg',files(ifile).folder);
        if ~exist(destination_folder), mkdir (destination_folder), end
       
        source = sprintf('%s/%s',files(ifile).folder,files(ifile).name);
        destination = sprintf('%s/%s',destination_folder,files(ifile).name);
        
        movefile(source,destination)
        
    end
    
    
    
end
