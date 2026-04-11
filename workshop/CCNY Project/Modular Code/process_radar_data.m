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
    
    % 4. DBSCAN Clustering
    [row, col] = find(detMapMasked > 0);
    if ~isempty(row)
        idx = dbscan([row, col], 5, 1); % min_pts=1 for high sensitivity
        numClusters = max(idx);
        targets = zeros(numClusters, 2); % Preallocate for [range, power]

        for k = 1:numClusters
            % Logic to find the cluster points
            clusterMask = (idx == k);
            clusterPtsR = row(clusterMask);
            
            % Calculate Centroid and Max Power
            centroidR = mean(clusterPtsR);
            clusterPower = max(magnitudeMTIMasked(clusterMask));
            
            % Assign directly to the preallocated row
            targets(k, :) = [round(centroidR), clusterPower];
        end
    end
end


