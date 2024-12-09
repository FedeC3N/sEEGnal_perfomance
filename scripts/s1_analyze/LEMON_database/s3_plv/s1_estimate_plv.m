clear
clc
restoredefaultpath

% Add functions to the path
addpath('../shared/')
addpath('../../../../../SharedFunctions/fieldtrip-20220626/')
ft_defaults

% Paths
config.path.clean_data = '../../../../databases/LEMON_database/derivatives';

% Define the frequency bands
bands_info = struct('name',{'delta', 'theta','alpha','beta','gamma'},...
    'f_limits',{[ 2 4] , [4 8] , [8 12] , [12 20] , [20 45] });

% To define later the PLV matrix
complete_channel_labels = {'Fp1', 'Fp2', 'F7', 'F3', 'Fz', 'F4', 'F8', 'FC5', 'FC1', 'FC2', 'FC6', 'T7', 'C3', 'Cz', 'C4', 'T8', 'CP5', 'CP1', 'CP2', 'CP6', 'AFz', 'P7', 'P3', 'Pz', 'P4', 'P8', 'PO9', 'O1', 'Oz', 'O2', 'PO10', 'AF7', 'AF3', 'AF4', 'AF8', 'F5', 'F1', 'F2', 'F6', 'FT7', 'FC3', 'FC4', 'FT8', 'C5', 'C1', 'C2', 'C6', 'TP7', 'CP3', 'CPz', 'CP4', 'TP8', 'P5', 'P1', 'P2', 'P6', 'PO7', 'PO3', 'POz', 'PO4', 'PO8'};

% Get the different testers
testers = dir(sprintf('%s/*',config.path.clean_data));
testers = testers(3:end);

for itester = 1 : numel(testers)
    
    % Load the datset
    current_tester = testers(itester);
    dataset_path = sprintf('%s/%s/%s_dataset.mat',current_tester.folder,...
        current_tester.name,current_tester.name);
    load(dataset_path);
    
    for ifile = 1 : numel(dataset)
        
        fprintf('Working on %s\n', dataset(ifile).file)
        
        % Process the subject
        [~,data] = evalc( 'process_subject(dataset(ifile))' );
        data_matrix = cat(3,data.trial{:});
        
        % Force to have 4 seconds (4096 samples)
        n = 1:4*data.hdr.Fs;
        data_matrix = data_matrix(:,n,:);
        
        % Create an index of the channels position
        channels_included = data.label;
        channels_included_index = ismember(complete_channel_labels,channels_included);
        
        % Output struct
        plv = [];
        plv.bands_info = bands_info;
        plv.channels_included = channels_included;
        plv.channels_included_index = channels_included_index;
        
        % Go through each band
        for iband = 1 : numel(bands_info)
            
            fprintf('  Band %s. ', bands_info(iband).name)
            
            % Filters the data in the selected band.
            band_order = 1800;
            current_f_limits = bands_info(iband).f_limits;
            current_banddata = myfiltfilt(data_matrix,band_order,current_f_limits,data.hdr.Fs);
            
            % Gets the data dimensions.
            nchans       = numel ( complete_channel_labels);
            nsamples     = size ( current_banddata, 2 );
            ntrials      = size ( current_banddata, 3 );
            
            fprintf ( 1, '  Calculating sensor-space connectivity.\n' );
            
            % Underlying logic:
            % The more time-consuming calculation is the complex exponential.
            % This algorithm works in the exponential space the whole time.
            % Uses the logarithmic operations for the phase difference:
            % exp ( 1i * ( a - b ) ) = exp ( 1i * a ) / exp ( 1i * b )
            % Performs all the divisions at once with matrix multiplication.
            
            % Memory reservation.
            PLV_all      = complex ( nan ( nchans, nchans, ntrials, 'single' ) );
            
            % Operates for each trial.
            for tindex = 1: ntrials
                
                %%% Calculate the complex PLV
                % Gets the sources time-series for the current trial.
                sensordata   = current_banddata ( :, :, tindex );
                
                % Gets the vector of phases for each time series.
                sensorvector = sensordata ./ abs ( sensordata );
                
                % Calculates the PLV as a matrix multiplication.
                PLV_all ( channels_included_index, channels_included_index, tindex )   = sensorvector * sensorvector' / nsamples;
                
                
                clear sourcedata sourcevector
                clear cPLV_t iplv_t rplv_t iplv_rms_t rplv_rms_t
            end
            
            % Sets the all-sources diagonal to NaN.
            PLV_all ( repmat ( eye ( nchans ) == 1, [ 1 1 ntrials ] ) ) = nan;
            
            % Average across trials
            PLV_all = mean(PLV_all,3);
            
            % Store the current PLV
            plv.(bands_info(iband).name).plv = PLV_all;
            
            clear PLV_all
        end
        
        % Save the PLV file
        path_plv = sprintf('../../../../databases/LEMON_database/derivatives/%s/plv',...
            current_tester.name);
        if ~exist(path_plv), mkdir(path_plv), end
        outfile_name = sprintf('sub-%s_ses-%s_task-%s_%s_plv',...
            dataset(ifile).sub,dataset(ifile).ses,dataset(ifile).task,...
            dataset(ifile).origin);
        outfile = sprintf('%s/%s',path_plv,outfile_name);
        save(outfile,'-struct','plv')
        
        % Save the metadata in the dataset
        plv = [];
        plv.path = path_plv;
        plv.file = outfile_name;
        dataset(ifile).plv = plv;
        
        fprintf('\n\n')
        
    end
    
    % Save
    save('-v7.3',dataset_path,'dataset')
    
    
end


% Functions
function  band_data = myfiltfilt(original_data,band_order,f_limits,fs);

% Creates the digital filter
b = fir1 ( band_order, f_limits / ( fs / 2 ) );

% Add padding
padding = size(original_data,2)/2;
band_data = padarray(original_data,[0 padding 0],'symmetric','both');

% Filter each channel and trial
for itrial = 1 : size(band_data,3)
    
    for ichannel = 1 : size(band_data,1)
        
        % Filter the data
        current_data = band_data(ichannel,:,itrial);
        current_data = filtfilt(b,1,current_data);
        
        % Put the data back
        band_data(ichannel,:,itrial) = current_data;
        
    end
    
end

% Removes the padding
band_data = band_data(:,padding+1:end-padding,:);

end
