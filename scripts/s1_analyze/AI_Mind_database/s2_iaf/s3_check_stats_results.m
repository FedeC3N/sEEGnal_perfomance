clear
clc
restoredefaultpath

% Paths
config.path.stats = '../../../../data/AI_Mind_database/stats';

% Load the results
load(sprintf('%s/iaf_stats.mat',config.path.stats));

% To export the result
txt_file = sprintf('%s/iaf_significant_results.txt',config.path.stats);
if exist(txt_file), delete(txt_file), end

% Table to excel
varNames = {'Measure','mean_eeg_expert(std_mean_err)','mean_ETL(std_mean_error)','t','p','Cohen d'};
varTypes = {'string','string','string','double','double','double'};
results_table = table('Size',[1,6],'VariableNames',varNames,'VariableTypes',varTypes);
counter = 1;

% Print the significant results for posterior analysis
measures = {'iaf','iaf_amp'};

% Stats
for imeasure = 1 : numel(measures)
    
    current_measure = measures{imeasure};
    current_mean_eeg_expert = sprintf('%.3f(%.3f)',stats.(current_measure).mean_eeg_expert,stats.(current_measure).std_mean_error_eeg_expert);
    current_mean_ETL = sprintf('%.3f(%.3f)',stats.(current_measure).mean_ETL,stats.(current_measure).std_mean_error_ETL);
    current_p = stats.(current_measure).p;
    current_t = stats.(current_measure).stats.tstat;
    current_d = stats.(current_measure).effect_size;
    
    % Add to table
    results_table(counter,:) = {current_measure,current_mean_eeg_expert,current_mean_ETL,...
        current_t,current_p, current_d};
    counter = counter + 1;
    
    % If significant, plot and save it in txt
    if current_p < 0.05
        
        
        line = sprintf('%s : p = %.2f\n', current_measure,current_p);
        
        % To screen
        fprintf(1,line);
        
        % To txt
        fid = fopen(txt_file,'a');
        fprintf(fid,line);
        fclose(fid);
        
    end
end




% Write the table to Excel
excel_file = sprintf('%s/iaf_significant_results.xlsx',config.path.stats);
if exist(excel_file), delete(excel_file),end
writetable(results_table,excel_file,'Sheet',1,'Range','A1')

