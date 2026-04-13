function [targets, detMapMasked, magnitudeMTIMasked] = process_radar_data(data, rdFiltered, mtiKernel, maskRange, cfar_thresh)
    % 1. MTI Filtering
    ncoeff = length(mtiKernel);
    mtiData = filter(mtiKernel, 1, data, [], 2);
    mtiData = mtiData(:, ncoeff:end);
    
    % 2. Range-Doppler Response
    [respMTI, ~, ~] = rdFiltered(mtiData);
    magnitudeMTI_Full = abs(respMTI);
    
    % 3. CFAR Detection
    [detMapFull, ~] = cfar_filter(magnitudeMTI_Full, 5, 5, 10, 8, cfar_thresh);
    detMapMasked = detMapFull(maskRange, :);
    magnitudeMTIMasked = magnitudeMTI_Full(maskRange, :);
    
    % 4. DBSCAN Clustering & Feature Extraction
    [row, col] = find(detMapMasked > 0);
    if ~isempty(row)
        idx = dbscan([row, col], 5, 1); % min_pts=1 for high sensitivity
        numClusters = max(idx);
        
        % PREALLOCATE FOR 7 FEATURES: 
        % [Average_Range,Doppler_Spread, Average_Doppler_Frequency, totalPower, Peak_to_Average_Power_Ratio, powerVariance, powerKurtosis]
        targets = zeros(numClusters, 7); 

        for k = 1:numClusters
            % 1. Isolate the cluster's specific points
            clusterMask = (idx == k);
            clusterPtsR = row(clusterMask);               % Range bins
            clusterPtsC = col(clusterMask);               % Doppler bins
            clusterPwr = magnitudeMTIMasked(clusterMask); % Power values
            
            % 2. Feature Extraction
            Average_Range = mean(clusterPtsR);
            totalPower = sum(clusterPwr);
            peakPower = max(clusterPwr);
            meanPower = mean(clusterPwr);
            Peak_to_Average_Power_Ratio = peakPower / meanPower;
            powerVariance = var(clusterPwr);
            Doppler_Spread = max(clusterPtsC) - min(clusterPtsC);

            powerKurtosis = kurtosis(clusterPwr);
            normPwr = clusterPwr / totalPower; % Normalize weights
            Average_Doppler_Frequency = sum(clusterPtsC .* normPwr); % Power-weighted centroid
            
            % 5. Assign the feature vector to the output row
            targets(k, :) = [Average_Range,Doppler_Spread, Average_Doppler_Frequency, totalPower, Peak_to_Average_Power_Ratio, powerVariance, powerKurtosis];
        end
    else
        targets = []; % Return empty if no targets pass CFAR
    end
end


