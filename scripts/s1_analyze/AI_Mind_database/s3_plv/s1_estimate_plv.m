clear
clc
restoredefaultpath

% Add functions to the path
addpath('../shared/')
addpath('../../../../../../SharedFunctions/fieldtrip-20220626/')
ft_defaults

% Paths
config.path.dataset = '../../../../data/SRM_database/dataset';
config.path.plv = '../../../../data/SRM_database/plv';

if ~exist(config.path.plv), mkdir(config.path.plv), end

% Define the frequency bands
% Define the frequency bands
bands_info = struct('name',{'delta', 'theta','alpha','beta','gamma'},...
    'f_limits',{[ 2 4] , [4 8] , [8 12] , [12 20] , [20 45] });

% To define later the PLV matrix
complete_channel_labels = {'Fp1';'AF7';'AF3';'F1';'F3';'F5';'F7';'FT7';'FC5';'FC3';'FC1';'C1';'C3';'C5';'T7';'TP7';'CP5';'CP3';'CP1';'P1';'P3';'P5';'P7';'P9';'PO7';'PO3';'O1';'Iz';'Oz';'POz';'Pz';'CPz';'Fpz';'Fp2';'AF8';'AF4';'AFz';'Fz';'F2';'F4';'F6';'F8';'FT8';'FC6';'FC4';'FC2';'FCz';'Cz';'C2';'C4';'C6';'T8';'TP8';'CP6';'CP4';'CP2';'P2';'P4';'P6';'P8';'P10';'PO8';'PO4';'O2'};

% Load the whole dataset
load(sprintf('%s/SRM_dataset.mat',config.path.dataset));

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
    outfile = sprintf('%s/sub-%s_ses-%s_task-%s_%s_plv', config.path.plv,...
        dataset(ifile).sub,dataset(ifile).ses,dataset(ifile).task,...
        dataset(ifile).database);
    save(outfile,'-struct','plv')
    
    % Save the metadata in the dataset
    plv = [];
    plv.path = config.path.plv;
    plv.file = sprintf('sub-%s_ses-%s_task-%s_%s_plv',...
        dataset(ifile).sub,dataset(ifile).ses,dataset(ifile).task,...
        dataset(ifile).database);
    dataset(ifile).plv = plv;
    
    fprintf('\n\n')
    
end

% Save
outfile = sprintf('%s/SRM_dataset.mat',config.path.dataset);
save('-v7.3',outfile,'dataset')



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
