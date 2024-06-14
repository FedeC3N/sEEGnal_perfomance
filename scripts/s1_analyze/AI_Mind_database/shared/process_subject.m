function data = process_subject(current_dataset)
    
    dummy_complete_file = sprintf('%s/%s',current_dataset.path,...
        current_dataset.file);
    
    % Load the file
    cfg = [];
    cfg.dataset = dummy_complete_file;
    data = ft_preprocessing(cfg);


end