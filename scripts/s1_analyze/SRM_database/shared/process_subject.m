function data = process_subject(current_dataset)
    
    dummy_complete_file = sprintf('%s/%s',current_dataset.path,...
        current_dataset.file);
    
    % Load the file
    cfg = [];
    cfg.dataset = dummy_complete_file;
    data = ft_preprocessing(cfg);
    
    % Process the data
    switch current_dataset.database
        
        case 'SRM_database'
            
            % Read the badchannels
            badchannel_file = dir(sprintf('%s/*channels.tsv',...
                current_dataset.path));
            badchannel_file = sprintf('%s/%s',current_dataset.path,...
                badchannel_file.name);
            
            % Fields of interest: name and status
            tsv_content = tdfread(badchannel_file);
            
            % Find the badchannels
            bad_index = ismember(cellstr(tsv_content.status),'bad');
            badchannels = cellstr(tsv_content.name);
            badchannels = badchannels(bad_index);
            
            % Remove the badchannels
            cfg = [];
            badchannels          = strcat('-',badchannels);
            desired_channel      = cat(2,{'eeg'},badchannels');
            cfg.channel          = desired_channel;
            cfg.precision        = 'single';
            cfg.feedback         = 'no';

            data            = ft_preprocessing ( cfg, data );
            
            
    end


end