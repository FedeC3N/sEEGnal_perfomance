clear
clc
restoredefaultpath

% Paths
config.path.stats = '../../../../data/AI_Mind_database/stats';

% Load the results
load(sprintf('%s/pow_stats.mat',config.path.stats));

% To export the result
txt_file = sprintf('%s/pow_significant_results.txt',config.path.stats);
if exist(txt_file), delete(txt_file), end

% Table to excel
varNames = {'Band','Area','mean_SRM (std_mean_err)','mean_ETL (std_mean_error)',...
    't','p','Cohen d'};
varTypes = {'string','string','string','string','double','double','double'};
results_table = table('Size',[1,7],'VariableNames',varNames,'VariableTypes',varTypes);
counter = 1;

% Print the significant results for posterior analysis
for iband = 1 : numel(bands_info)
   
    current_band = bands_info(iband).name;
    
    for iarea = 1 : numel(areas_info)
    
    	current_area = areas_info(iarea).name;
        
        % Stats
        current_mean_eeg_expert = sprintf('%10e (%10e)',stats.(current_band).(current_area).mean_eeg_expert,stats.(current_band).(current_area).std_mean_error_eeg_expert);
        current_mean_ETL = sprintf('%10e (%10ef)',stats.(current_band).(current_area).mean_ETL,stats.(current_band).(current_area).std_mean_error_ETL); 
        current_p = stats.(current_band).(current_area).p;
        current_t = stats.(current_band).(current_area).stats.tstat;
        current_d = stats.(current_band).(current_area).effect_size;
            
        % Add to table
        if iarea == 1
            results_table(counter,:) = {current_band, current_area,...
                current_mean_eeg_expert,current_mean_ETL,current_t,current_p, current_d};
        else
            results_table(counter,:) = {'', current_area,...
                current_mean_eeg_expert,current_mean_ETL,current_t,current_p, current_d};
        end
        counter = counter + 1;
        
        % If significant, plot and save it in txt
        if current_p < 0.05
            
            
            line = sprintf('%s - %s : p = %.2f\n', current_band, current_area,current_p);
            
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
excel_file = sprintf('%s/pow_significant_results.xlsx',config.path.stats);
if exist(excel_file), delete(excel_file),end
writetable(results_table,excel_file,'Sheet',1,'Range','A1')

