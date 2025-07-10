%{

Wrapper to use Fieldtrip to read the EEG files.
Re-reference the recordings to the average.

@author: Fede

%}

function data = process_subject(current_dataset)
    
    dummy_complete_file = fullfile('..','..','..','..',...
        current_dataset.path, current_dataset.file);
    
    cfg = [];
    cfg.dataset = dummy_complete_file;
    cfg.reref = 'yes';
    cfg.refmethod = 'avg';
    cfg.refchannel = 'all';
    data = ft_preprocessing(cfg);
    
    
    % Fix the header
%     hdr = ft_read_header(dummy_complete_file);
%     hdr.orig.filename = current_dataset.file;
%     hdr.orig.data = current_dataset.file;
%     hdr.orig.data(end-3:end) = '.set';
%     hdr.orig.datfile = hdr.orig.data;
%     data = ft_read_data(dummy_complete_file,'header',hdr);

end