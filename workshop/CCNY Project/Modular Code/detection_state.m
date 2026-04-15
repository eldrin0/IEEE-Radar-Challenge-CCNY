function [color, txt, stateOut] = detection_state(targets, magRaw, rGrid, currentTime)
    % Persistent variables keep their value between function calls
    persistent lastDetectionTime isTargetLocked lockedRangeIdx storedColor storedText
    
    % Initialization on first run
    if isempty(lastDetectionTime)
        lastDetectionTime = -inf;
        isTargetLocked = false;
        lockedRangeIdx = 1;
        storedColor = [0.8 0.8 0.8];
        storedText = 'Searching...';
    end

    % Parameters
    persistenceDuration = 1; % Time to wait before going back to gray
    rawThresholdFactor = 1;  % Stationary sensitivity (IMPORTANT)

    if ~isempty(targets)
        % --- PATH A: MOTION ACQUISITION ---
        % targets = [RangeIdx, Power]

        lockedRangeIdx = targets(1, 1);
        realRange = rGrid(lockedRangeIdx);
        lastDetectionTime = currentTime;
        isTargetLocked = true;
        
        if realRange <= 1.5
            storedColor = [1 0 0]; % Red
            storedText = sprintf('DANGER: %.2fm', realRange);
        elseif realRange <= 3.5
            storedColor = [1 0.9 0]; % Yellow
            storedText = sprintf('WARNING: %.2fm', realRange);
        else
            storedColor = [0 0.8 0]; % Green
            storedText = sprintf('CLEAR: %.2fm', realRange);
        end
        
    elseif isTargetLocked && (currentTime - lastDetectionTime) < 8.0
        % --- PATH B: STATIONARY HOLD ---
        % Look at the 0-Doppler energy at the last known range
        zeroDopplerIdx = floor(size(magRaw, 2) / 2) + 1;
        currentRawNoiseFloor = mean(magRaw(:));
        currentRawEnergy = magRaw(lockedRangeIdx, zeroDopplerIdx);
        
        if currentRawEnergy > (currentRawNoiseFloor * rawThresholdFactor)
            if ~contains(storedText, '(STATIONARY)')
                storedText = [storedText, ' (STATIONARY)'];
            end
            lastDetectionTime = currentTime; 
        else
            isTargetLocked = false;
        end
    end
    
    % --- TIMEOUT CHECK ---
    if (currentTime - lastDetectionTime) > persistenceDuration
        storedColor = [0.8 0.8 0.8]; 
        storedText = 'Searching...';
        isTargetLocked = false;
    end

    % Return values for the UI
    color = storedColor;
    txt = storedText;
    stateOut.isLocked = isTargetLocked; % Optional: for debugging
end