clear
clc
restoredefaultpath

% Paths
config.path.results = '../../../../results/plv';

% Load the results
load(sprintf('%s/plv_results.mat',config.path.results));

% To export the result
txt_file = sprintf('%s/plv_significant_results.txt',config.path.results);
if exist(txt_file), delete(txt_file), end

% Table to excel
varNames = {'Band','Area','mean_lemon (std_mean_err)','mean_sEEGnal (std_mean_error)',...
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
        current_mean_lemon = sprintf('%10e (%10e)',stats.(current_band).(current_area).mean_lemon,stats.(current_band).(current_area).std_mean_error_lemon);
        current_mean_sEEGnal = sprintf('%10e (%10ef)',stats.(current_band).(current_area).mean_sEEGnal,stats.(current_band).(current_area).std_mean_error_sEEGnal); 
        current_p = stats.(current_band).(current_area).p;
        current_t = stats.(current_band).(current_area).stats.tstat;
        current_d = stats.(current_band).(current_area).effect_size;
            
        % Add to table
        if iarea == 1
            results_table(counter,:) = {current_band, current_area,...
                current_mean_lemon,current_mean_sEEGnal,current_t,current_p, current_d};
        else
            results_table(counter,:) = {'', current_area,...
                current_mean_lemon,current_mean_sEEGnal,current_t,current_p, current_d};
        end
        counter = counter + 1;
        
        % If significant, plot and save it in txt
        if current_p < 0.05
            
            
            line = sprintf('%s - %s : p = %.5f\n', current_band, current_area,current_p);
            
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


