clear
clc
restoredefaultpath

% Paths
config.path.results = '../../../../results/LEMON_database/plv';

% Load the results
load(sprintf('%s/plv_results.mat',config.path.results));

% To export the result
txt_file = sprintf('%s/plv_significant_results.txt',config.path.results);
if exist(txt_file), delete(txt_file), end

% Table to excel
varNames = {'Band','Channel','t_mean','t_std','d_mean','d_std','p_mean','p_std'};
varTypes = {'string','string','double', 'double','double', 'double','double', 'double'};
results_table = table('Size',[1,8],'VariableNames',varNames,'VariableTypes',varTypes);
counter = 1;

% Print the significant results for posterior analysis
for iband = 1 : numel(bands_info)
   
    current_band = bands_info(iband).name;
    
    for ichannel = 1 : numel(complete_channel_labels)
    
        current_channel = complete_channel_labels{ichannel};
        
        % Stats
        current_t_mean = nanmean(stats.(current_band).tstat(ichannel,:));
        current_t_std = nanstd(stats.(current_band).tstat(ichannel,:));
        current_d_mean = nanmean(stats.(current_band).cohen_d(ichannel,:));
        current_d_std = nanstd(stats.(current_band).cohen_d(ichannel,:));
        current_p_mean = nanmean(stats.(current_band).p(ichannel,:));
        current_p_std = nanstd(stats.(current_band).p(ichannel,:));
            
        % Add to table
        if ichannel == 1
            results_table(counter,:) = {current_band, current_channel,...
                current_t_mean,current_t_std,current_d_mean,current_d_std,...
                current_p_mean,current_p_std};
        else
            results_table(counter,:) = {current_band, current_channel,...
                current_t_mean,current_t_std,current_d_mean,current_d_std,...
                current_p_mean,current_p_std};
        end
        counter = counter + 1;
        
        % If significant, plot and save it in txt
        if current_p_mean < 0.05
            
            
            line = sprintf('%s - %s : p = %.2f\n', current_band, current_channel,current_p_mean);
            
            % To screen
            fprintf(1,line);
            
            % To txt
            fid = fopen(txt_file,'a');
            fprintf(fid,line);
            fclose(fid);
            
        end
        
    end
    
    fprintf('\n')
    
    
end

% Write the table to Excel
excel_file = sprintf('%s/plv_significant_results.xlsx',config.path.results);
if exist(excel_file), delete(excel_file),end
writetable(results_table,excel_file,'Sheet',1,'Range','A1')

