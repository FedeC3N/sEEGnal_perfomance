function [pow_dataset,f] = read_pow_dataset(dataset, desired_dataset)

% subject of interest
current_dataset_index = ismember({dataset.database},desired_dataset);
current_dataset = dataset(current_dataset_index);

for icurrent = 1 : numel(current_dataset)
    
    % Load pow
    pow = load(sprintf('%s/%s',current_dataset(icurrent).pow.path,...
        current_dataset(icurrent).pow.file));
    
    % Save f
    f = pow.f;
    
    % Normalize
    current_pow = mean(pow.pow_spectrum,3);
    scaling_factor = 1./sum(current_pow,2);
    scaling_factor = repmat(scaling_factor,1,size(current_pow,2));
    current_pow_norm = scaling_factor .* current_pow;
    
    % Add to the all matrix
    if icurrent == 1
        pow_dataset = nan(64,numel(f),numel(current_dataset));
    end
    pow_dataset(1:size(current_pow_norm,1),:,icurrent) = current_pow_norm;
    
end

end