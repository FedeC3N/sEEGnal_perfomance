% Estimate Individual Alpha Frequency (and amplitude) in occipital sensors
clear
clc
restoredefaultpath

% Add functions to the path
addpath('../shared/')

% Paths
config.path.dataset = '../../../../data/SRM_database/dataset';
config.path.iaf = '../../../../data/SRM_database/iaf';

if ~exist(config.path.iaf), mkdir(config.path.iaf), end

% To define later the pow matrix
complete_channel_labels = {'Fp1';'AF7';'AF3';'F1';'F3';'F5';'F7';'FT7';'FC5';'FC3';'FC1';'C1';'C3';'C5';'T7';'TP7';'CP5';'CP3';'CP1';'P1';'P3';'P5';'P7';'P9';'PO7';'PO3';'O1';'Iz';'Oz';'POz';'Pz';'CPz';'Fpz';'Fp2';'AF8';'AF4';'AFz';'Fz';'F2';'F4';'F6';'F8';'FT8';'FC6';'FC4';'FC2';'FCz';'Cz';'C2';'C4';'C6';'T8';'TP8';'CP6';'CP4';'CP2';'P2';'P4';'P6';'P8';'P10';'PO8';'PO4';'O2'};
occipital_channels = {'O1';'Iz';'Oz';'O9';'O2';'O10'};

% Load the whole dataset
load(sprintf('%s/SRM_dataset.mat',config.path.dataset));


for ifile = 1 : numel(dataset)
    
    fprintf('Working on %s\n\n', dataset(ifile).file)
    
    % Load pow
    pow = load(sprintf('%s/%s',dataset(ifile).pow.path,...
        dataset(ifile).pow.file));
    
    % Channels of interest
    desired_channels_index = ismember(complete_channel_labels, occipital_channels);
    
    % Get the pow norm
    % Normalize
    current_pow = mean(pow.pow_spectrum,3);
    scaling_factor = 1./nansum(nansum(current_pow,2));
    scaling_factor = repmat(scaling_factor,size(current_pow));
    current_pow_norm = scaling_factor .* current_pow;
    
    % Get the channels of interest
    current_pow_norm = current_pow_norm(desired_channels_index,:);
    
    % Average the occipital channels
    current_pow_norm = nanmean(current_pow_norm,1);
    
    % If all occipital channels are badchannels, is a nan
    if isnan(current_pow_norm(1))
        peak = [];
        peak.iaf = nan;
        peak.iaf_amp = nan;
        
    else
        
        % Search the alpha peak
        f = pow.f';
        foi = f(f > 4 & f < 15)';
        foi = double(foi);
        options = optimset('Display','off');
        fitrange = [4 15]; % search_1peak would search within these frequencies.
        peak = search_alpha_peak(double(f),double(current_pow_norm),...
            double(fitrange),options);
        
        log_fitted_powspctrm = peak.B - ...
            peak.C*log(f) + ...
            peak.A*exp(-((f-peak.iaf)/peak.sp).^2);
        
        % Plot
        %     plot(f, log(current_pow_norm),'r'); hold on;
        %     plot(f, log_fitted_powspctrm, 'b');
        %     plot(peak.iaf, peak.iaf_amp,'k*')
        %     legend('Real log(powsprtrm)', 'Fitted powsprtrm' , 'Frequency peak');
        
    end
    
    % Add to the dataset struct
    iaf = [];
    iaf.peak_estimation = peak;
    iaf.peak_estimation.formula = 'log(pow)= B - C*log(f) + A*exp(-((f-fp)/sp).^2)';
    iaf.peak_estimation.f = f;
    iaf.iaf = peak.iaf;
    iaf.iaf_amp = peak.iaf_amp;
    
    % Save the file
    outfile = sprintf('%s/sub-%s_ses-%s_task-%s_%s_iaf', config.path.iaf,...
        dataset(ifile).sub,dataset(ifile).ses,dataset(ifile).task,...
        dataset(ifile).database);
    save(outfile,'-struct','iaf')
    
    % Save the metadata in the dataset
    iaf = [];
    iaf.path = config.path.iaf;
    iaf.file = sprintf('sub-%s_ses-%s_task-%s_%s_iaf',...
        dataset(ifile).sub,dataset(ifile).ses,dataset(ifile).task,...
        dataset(ifile).database);
    dataset(ifile).iaf = iaf;
    
end

% Save
outfile = sprintf('%s/SRM_dataset.mat',config.path.dataset);
save('-v7.3',outfile,'dataset')


% Functions
function peak = search_alpha_peak(foi,data,rangepeak,options)

% function that fits a peak with power law background to the powerspectrum
% fits according to log(pow)= B - C*log(f) + A*exp(-((f-fp)/sp).^2)
% and finds the corresponding parameters B, C, A, fp and sp

% INPUTS
%   foi= 1x nfoi frequencies of interest
%   data= ntrials x nfoi = powerspectrum for each foi and trial. The peak
%      can also be fitted to a single or average trial, but in that case
%      the intertrial variability cannot be estimated and the output
%      peakfittedok is meaningless
%   rangepeak = [fmin fmax] a-priori range (in Hz) for the peak (ex: [5-14])


% OUTPUTS
%   peaks: structure containing the values of the fitted parameters B, C,
%       A, fp and sp
%   peakfittedok: true or false, depending on whether the amplitude of the
%       peak exceeds the rms intertrial variability


% implemented according to recommendations in:
%    Lodder, S. S., & van Putten, M. J. (2011). Automated EEG analysis:
%        Characterizing the posterior dominant rhythm. Journal of
%        neuroscience methods, 200(1), 86-93.
%    Chiang, A. K. I., Rennie, C. J., Robinson, P. A., Roberts, J. A.,
%        Rigozzi, M. K., Whitehouse, R. W., ... & Gordon, E. (2008).
%        Automated characterization of multiple alpha peaks in multi-site
%        electroencephalograms. Journal of neuroscience methods, 168(2),
%        396-411.


peak=[];
logpf=log(trimmean(data,10,'round',1)); %Average over trials --> logpf: 1 x ntrials
Ntrials=size(data,1);

idrangepeak=find(and(foi>rangepeak(1),foi<rangepeak(2))); %range for the peak frequency


% 1 - Estimates background spectrum logpbackground = B - C*log(f)
pfminusbackground = @(p)(logpf + p(1)*log(foi) - p(2)  );
minB=polyfit(foi,logpf,0);
maxB=polyfit(foi,logpf+2*log(foi),0);
[p,~] = lsqnonlin(pfminusbackground,double([0 mean(logpf)]),...
    double([0 minB]),double([2 maxB]),options);
peak.C=p(1); peak.B=p(2);

% 2 - Fits peak

% 2a-initial value for peak frequency and amplitude in the optimization algorithm
[vmax,imax]=max(logpf(idrangepeak)-(peak.B - peak.C*log(foi(idrangepeak))));
f1=foi(idrangepeak(imax));

% 2b - finds peak parameters with lsqnonlin optimization
% C stays fixed from the previous optimization
myerror = @(x)(abs(x(1) - peak.C*log(foi) + x(2)*exp(-((foi-x(4))/x(3)).^2) - logpf));
[x,~] = lsqnonlin(myerror,[peak.B; vmax; 1; f1],[minB; 0; 0.5; rangepeak(1)],[maxB; 100*vmax; 3; rangepeak(2)],options);
peak.B=x(1); peak.A=x(2); peak.sp=x(3); peak.iaf=x(4);
peak.iaf_amp = peak.B - peak.C*log(peak.iaf) + peak.A;



end