clear
clc

path_badchannels = '../../../aimind.etl.data/ETL_performance_AIMIND_dataset/standardized/eeg/derivatives/etl';
path_out = '../../../aimind.etl.data/ETL_performance_AIMIND_dataset/metadata/unified_format';

subjects = dir(sprintf('%s/*',path_badchannels));
subjects = {subjects.name};
subjects = subjects(3:end);

recording_type = {'EO' 'EC'};

load('channels_labels.mat');

for isubject = 1:numel(subjects)
    
    current_subject = subjects{isubject};
    
    for itype = 1:numel(recording_type)
        
        current_type = recording_type{itype};
        
        % Bad channels
        files = dir(sprintf('%s/%s/ses-1/eeg/%s*%s*channels.tsv',...
            path_badchannels,current_subject,current_subject,current_type));
        current_file = sprintf('%s/%s/ses-1/eeg/%s',...
            path_badchannels,current_subject,files.name);
        
        s = tdfread(current_file);
        
        if isfield(s,'bad')
            bads = cat(1,'bad ',s.bad);
        else
            bads = cat(1,'good',s.good);
        end
        bads = cellstr(bads);
        bads_mask = strcmp(bads,'bad');
        
        % Output structure
        output = [];
        output.badchannels.channel_labels = channel_labels;
        output.badchannels.badchannels_labels = channel_labels(bads_mask);
        output.badchannels.badchannels_index = find(bads_mask);
        
        dummy = files.name;
        dummy = dummy(1:end-13);
        outfile = sprintf('%s/IClabel_processed_%s',path_out,dummy);
        
        % annotations
        files = dir(sprintf('%s/%s/ses-1/eeg/%s*%s*_artifacts_annotations.tsv',...
            path_badchannels,current_subject,current_subject,current_type));
        files = {files.name};
        
        % Discard SOBI files
        index_file = ~contains(files,'sobi');
        files = files(index_file);
        
        current_file = sprintf('%s/%s/ses-1/eeg/%s',...
            path_badchannels,current_subject,files{1});
        
        s = tdfread(current_file);
        
        % Output structure
        output.annotations.artifact_type = cellstr(s.label);
        output.annotations.start_sample = s.onset * 2000;
        output.annotations.end_sample = (s.onset + s.duration)*2000;
          
        save(outfile,'-struct','output')
        
        
    end
        
        
end


