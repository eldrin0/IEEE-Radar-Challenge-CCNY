function heatmap = HeatMap_Collection(validRanges, nAngles, bf, precomputedPhases, rx, bf_TDD, digitalsteer, rr, keepRange)
    heatmapMatrix = zeros(length(validRanges), nAngles);
    
    % --- Apply Range Compensation (Time Variable Gain) ---
    % We add 0.1 to prevent multiplying Range=0 by absolute zero.
    % We use an exponent of 2 to compensate for the two-way voltage loss.
    range_weights = (validRanges(:) + 0.1).^2; 
    
    for i = 1:nAngles
        bf.RxPhase = precomputedPhases(i, :);
        bf.LatchRxSettings();
        
        bf.Burst = true; 
        data = rx();
        bf.Burst = false;
        
        dataProcessed = arrangePulseData(data, rx, bf, bf_TDD);
        
        if size(dataProcessed, 2) >= 2
            ch1 = dataProcessed(:, 1); 
            ch2 = dataProcessed(:, 2);
        else
            ch1 = dataProcessed(:, 1);
            ch2 = zeros(size(ch1)); 
        end
        
        dataCombined = (ch1 .* conj(digitalsteer(1))) + (ch2 .* conj(digitalsteer(2)));
        [resp, ~] = rr(dataCombined);
        
        % Extract magnitude and multiply by the range weights
        raw_mag = abs(resp(keepRange));
        heatmapMatrix(:, i) = raw_mag .* range_weights;
    end
    
    heatmap = heatmapMatrix;
end