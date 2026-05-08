function [lockedRangeIdx, lastDetectionTime, isTargetLocked, storedColor, storedText] = ...
    update_state(targets, lockedRangeIdx, rangeGridFiltered, currentTime, ...
                 lastDetectionTime, isTargetLocked, sGridRaw, magnitudeRaw, ...
                 storedColor, storedText)
             
    if ~isempty(targets)
        lockedRangeIdx = targets(1, 1);
        realRange = rangeGridFiltered(lockedRangeIdx);
        if realRange <= 1, currentZone = 3;
        elseif realRange <= 3, currentZone = 2;
        else, currentZone = 1; 
        end
        
        lastDetectionTime = currentTime;
        isTargetLocked = true;
        
        if currentZone == 3
            storedColor = [1 0 0]; storedText = sprintf('Stop: %.2fm', realRange);
        elseif currentZone == 2
            storedColor = [1 0.9 0]; storedText = sprintf('Slow Down: %.2fm', realRange);
        elseif currentZone == 1
            storedColor = [0 0.8 0]; storedText = sprintf('Warn: %.2fm', realRange);
        end
    elseif isTargetLocked && (currentTime - lastDetectionTime) < 5.0
        [~, zeroDopplerIdx] = min(abs(sGridRaw)); 
        currentRawNoiseFloor = mean(magnitudeRaw(:));
        currentRawEnergy = magnitudeRaw(lockedRangeIdx, zeroDopplerIdx);
        if currentRawEnergy > (currentRawNoiseFloor * 15)
            if ~contains(storedText, '(STATIONARY)'), storedText = [storedText, ' (STATIONARY)']; end
            lastDetectionTime = currentTime; 
        else
            isTargetLocked = false;
            storedColor = [0.8 0.8 0.8]; storedText = 'Clear';
        end
    end
    
    if (currentTime - lastDetectionTime) > 2.5
        storedColor = [0.8 0.8 0.8]; storedText = 'Clear';
        isTargetLocked = false;
    end
end