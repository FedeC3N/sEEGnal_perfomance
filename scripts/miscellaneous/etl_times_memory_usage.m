clear
clc
restoredefaultpath

path.metadata = '../../metadata/AI_Mind_database/trl/etl';

recording_type = {'EO', 'EC'};

for itype =  1 : numel(recording_type)
    
    files = dir(sprintf('%s/*%s*.mat',path.metadata,recording_type{itype}));
    files = {files.name};
    files = files(1:end-1);
    
    process = {'standardization','badchannel_detection','artifact_detection','final_qa'}';
    times = nan(numel(process),1);
    memory_MB = nan(numel(process),1);
    
    for iprocess = 1 : numel (process)
        
        current_times = nan(numel(files),1);
        current_memory_MB = nan(numel(files),1);
        
        for ifile = 1 : numel(files)
            
            trl = load(sprintf('%s/%s',path.metadata,files{ifile}));
            
            current_times(ifile) = trl.times_seconds.(process{iprocess});
            current_memory_MB(ifile) = trl.memory_bytes.(process{iprocess})/1024/1024;
            
        end
        
        times(iprocess) = nanmean(current_times);
        memory_MB(iprocess) = nanmean(current_memory_MB);
    end
    
    ETL_performance_table = table(process,times,memory_MB);
    fprintf(1,'\n %s recordings \n\n',recording_type{itype})
    disp(ETL_performance_table)
    fprintf(1,'\n\n')
end
