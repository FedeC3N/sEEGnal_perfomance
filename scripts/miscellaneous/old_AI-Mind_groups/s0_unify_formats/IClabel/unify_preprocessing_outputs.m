clear
clc

path_derivatives = '../../../data/AI-Mind_database/standardized/derivatives/etl';
path_out = '../../../data/AI-Mind_database/metadata/unified_format';

subjects = dir(sprintf('%s/sub*',path_derivatives));
subjects = {subjects.name};

recording_type = {'1EO' '2EC', '3EO' '4EC'};

load('../channels_label.mat');

%estructura_como_guardar = load(sprintf('%s/Fede_processed_1-103-1-Z_1-EO.mat',path_out));

for isubject = 1:numel(subjects)
    
    current_subject = subjects{isubject};
    
    sessions = dir(sprintf('%s/%s/ses*',path_derivatives,current_subject));
    sessions = {sessions.name};
    
    for isession = 1:numel(sessions)
        
        current_session = sessions{isession};
        
        
        for itype = 1:numel(recording_type)
            
            current_type = recording_type{itype};
            
            % Bad channels
            files = dir(sprintf('%s/%s/%s/eeg/%s*%s*channels.tsv',...
                path_derivatives,current_subject,current_session,...
                current_subject,current_type));
            if isempty(files)
                fprintf('No file. Check %s\n\n',current_subject);
                continue
            end
            current_file = sprintf('%s/%s/%s/eeg/%s',...
                path_derivatives,current_subject,current_session,files.name);
            
            s = tdfread(current_file);
                        
            bads = cellstr(s.status);
            bads_mask = strcmp(bads,'bad');
            
            channels_in_file = cellstr(s.x0xEF0xBB0xBFname);
            current_badchannel_labels = channels_in_file(bads_mask);
            current_badchannels_index = find(ismember(channels_label,current_badchannel_labels));
            
            % Output structure
            output = [];
            output.badchannels.channel_labels = channels_label;
            output.badchannels.badchannels_labels = current_badchannel_labels;
            output.badchannels.badchannels_index = current_badchannels_index;

            % annotations
            files = dir(sprintf('%s/%s/%s/eeg/%s*%s*artifacts_annotations.tsv',...
                path_derivatives,current_subject,current_session,...
                current_subject,current_type));
            files = {files.name};
            if isempty(files)
                fprintf('No file. Check %s\n\n',current_subject);
                continue
            end
            
            % Discard SOBI files
            index_file = ~contains(files,'sobi');
            files = files(index_file);
            
            current_file = sprintf('%s/%s/%s/eeg/%s',...
                path_derivatives,current_subject,current_session,files{1});
            
            try
                s = tdfread(current_file);
                
                % Output structure
                output.annotations.artifact_type = cellstr(s.label);
                output.annotations.start_sample = s.x0xEF0xBB0xBFonset * 2000;
                output.annotations.end_sample = (s.x0xEF0xBB0xBFonset + s.duration)*2000;
                
            catch
                
                 % Output structure
                output.annotations.artifact_type = [];
                output.annotations.start_sample = [];
                output.annotations.end_sample = [];
                
            end
            
            
            
            % IC
            % Read mixing matrix
            files = dir(sprintf('%s/%s/%s/eeg/%s*%s*sobi_mixing.tsv',...
                path_derivatives,current_subject,current_session,...
                current_subject,current_type));
            if isempty(files)
                fprintf('No file. Check %s\n\n',current_subject);
                continue
            end
            current_file = sprintf('%s/%s/%s/eeg/%s',...
                path_derivatives,current_subject,current_session,files.name);
            s = tdfread(current_file);
            
            % Construct the mixing matrix
            mixing = zeros(size(s.IC001,1));
            field_names = fieldnames(s);
            index_field_names = (regexp(field_names,'IC*','ONCE'));
            index_field_names = ~cellfun(@isempty, index_field_names);
            field_names = field_names(index_field_names);
            for ifield = 1:numel(field_names)
               
                current_field = field_names{ifield};               
                mixing (:,ifield) = s.(current_field);
                
            end
        
            % Construct the unmixing matrix
            unmixing = pinv(mixing);
            
            % excluded_mask
            files = dir(sprintf('%s/%s/%s/eeg/%s*%s*sobi_annotations.tsv',...
                path_derivatives,current_subject,current_session,...
                current_subject,current_type));
            if isempty(files)
                fprintf('No file. Check %s\n\n',current_subject);
                continue
            end
            current_file = sprintf('%s/%s/%s/eeg/%s',...
                path_derivatives,current_subject,current_session,files.name);
            s = tdfread(current_file);
            
            % Select brain + other
            excluded_mask = false(size(mixing,1),1);
            for ilabel = 1:size(s.label,1)
               
                current_label = s.label(ilabel,:);
                current_label = strtrim(current_label);
                
                if (strcmp('brain',current_label) || strcmp('other',current_label))
                    continue
                end
                
                current_excluded = s.channel(ilabel,:);
                current_excluded = strtrim(current_excluded);
                if strcmp(current_excluded,'n/a')
                    continue
                end
                current_excluded = split(current_excluded,',');
                current_excluded = cellfun(@(x) x(3:end),current_excluded,'UniformOutput',false);
                current_excluded = cellfun(@str2double,current_excluded,'UniformOutput',false);
                current_excluded = cell2mat(current_excluded);
                
                
                excluded_mask(current_excluded) = true;
            
            end
            IC_included_mask = ~excluded_mask;
            
            
            % Create the output IC
            ICs = [];
            ICs.mixing_matrix = mixing;
            ICs.unmixing_matrix = unmixing;
            ICs.IC_included_mask = IC_included_mask;
            
            % Save the dta
            output.ICs = ICs;
            
            % Create the file name
            dummy = current_session(5:end);
            dummy = sprintf('%s-%s-%s-%s_%s-%s',dummy(1),dummy(2:4),dummy(5),...
                dummy(6),current_type(1),current_type(2:3));

            outfile = sprintf('%s/IClabel_processed_%s',path_out,dummy);
            
            save(outfile,'-struct','output')
            
            
        end
        
    end
    
end


