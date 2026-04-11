clear all; close all; clc;
clear classes;

% --- 1. Parameters (Sensitivity Tuning) ---
cfar_threshold = 1.0e-2; % LOWERED: Was 3e-2. Lower = more sensitive.
min_points_cluster = 1;  % LOWERED: Was 3. 1 means "anything" is a target.
confirmThreshold = 1;    % LOWERED: Was 3. 1 means immediate UI update.

fc = 10e9; fs = 1e6; 
rampbandwidth = 2000e6; 
prf = 2000; tSweep = 5e-4; 
sweepSlope = rampbandwidth / tSweep;
calRange = 2.75; 
minRange = 0; maxRange = 10;
tCapture = 1200;
nPulses = 512;
mtiKernel = [1 -1]; 
ncoeff = length(mtiKernel);

% --- 2. Logic & Stability Initialization ---
confirmCount = 0;       
lastZone = 0;           
isTargetLocked = false;
lockedRangeIdx = 1;
storedColor = [0.8 0.8 0.8];
storedText = 'Searching...';
lastDetectionTime = -inf;

% --- 3. Hardware & Scope Setup ---
[rx,tx,bf,bf_TDD] = setupLabRadar(fc,prf,nPulses,fs,rampbandwidth);

rdRaw = phased.RangeDopplerResponse(DopplerOutput="Speed",...
    OperatingFrequency=fc,SampleRate=fs,RangeMethod="FFT",...
    SweepSlope=sweepSlope,PRFSource="Property",PRF=prf);

rdFiltered = phased.RangeDopplerResponse('DopplerOutput', 'Speed', 'OperatingFrequency', fc, ...
    'SampleRate', fs, 'RangeMethod', 'FFT', 'SweepSlope', sweepSlope, 'PRFSource', 'Property', 'PRF', prf);

scopeRaw = phased.RangeDopplerScope('IQDataInput', false, 'Name', '1. Raw Response', 'Position', [10 450 480 400]);
scopeMTI = phased.RangeDopplerScope('IQDataInput', false, 'Name', '2. MTI Filtered', 'Position', [500 450 480 400]);
scopeCFAR = phased.RangeDopplerScope('IQDataInput', false, 'Name', '3. CFAR Detections', 'Position', [990 450 480 400]);

% --- 3 & 4. Unified UI Setup ---
hMainFig = figure('Name', 'Radar Dashboard - Logic Fusion', 'Position', [100 100 1200 800], 'Color', 'w');
tlo = tiledlayout(2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

% Tile 1: Raw Response
axRaw = nexttile;
title(axRaw, '1. Raw Response');

% Tile 2: MTI Filtered
axMTI = nexttile;
title(axMTI, '2. MTI Filtered');

% Tile 3: CFAR Detections
axCFAR = nexttile;
title(axCFAR, '3. CFAR Detections');

% Tile 4: Diagnostic Alert System
axTarget = nexttile;
hold(axTarget, 'on');
hSquare = rectangle('Parent', axTarget, 'Position', [0.1, 0.2, 0.8, 0.6], ...
    'Curvature', [0.1, 0.1], 'FaceColor', [0.8 0.8 0.8], 'LineWidth', 2);
hText = text(axTarget, 0.5, 0.1, 'Searching...', 'HorizontalAlignment', 'center', ...
    'FontSize', 12, 'FontWeight', 'bold');
title(axTarget, 'Target Proximity Status');
set(axTarget, 'XTick', [], 'YTick', [], 'XColor', 'none', 'YColor', 'none');
axis(axTarget, [0 1 0 1]);

% --- Re-initialize Scopes to target the specific axes ---
% Note: Standard phased.RangeDopplerScope creates its own window. 
% To keep it in the tiles, we use standard imagesc or surf later in the loop.

% --- 5. Main Processing Loop ---
amp = 0.9 * 2^15;
txWaveform = amp * ones(rx.SamplesPerFrame, 2);
tx(txWaveform);
t_start = tic;

while toc(t_start) < tCapture
    currentTime = toc(t_start);
    
    % A. Data Capture
    bf.Burst = false; bf.Burst = true; bf.Burst = false;
    rawdata = rx();
    data = arrangePulseData(rawdata, rx, bf, bf_TDD);
   
    % B. Raw Processing
    [respRaw, rGridRaw, sGridRaw] = rdRaw(data);
    rGridRawCal = rGridRaw - calRange;
    maskRange = rGridRawCal >= minRange & rGridRawCal <= maxRange;
    magnitudeRaw = abs(respRaw(maskRange, :));
    rangeGridFiltered = rGridRawCal(maskRange);
    
    % C. MTI & CFAR Processing
    mtiData = filter(mtiKernel, 1, data, [], 2);
    mtiData = mtiData(:, ncoeff:end);
    [respMTI, rGridMTI, sGridMTI] = rdFiltered(mtiData);
    
    magnitudeMTI_Full = abs(respMTI);
    % Using the dynamic cfar_threshold variable from the top
    [detMapFull, noiseMapFull] = cfar_filter(magnitudeMTI_Full, 5, 5, 10, 8, cfar_threshold);
    
    detMapMasked = detMapFull(maskRange, :);
    magnitudeMTIMasked = magnitudeMTI_Full(maskRange, :);
    
    % Update Visuals
    % --- Update Visuals in Unified Window ---
    % Plot 1: Raw
    imagesc(axRaw, sGridRaw, rangeGridFiltered, 10*log10(magnitudeRaw));
    axis(axRaw, 'xy'); xlabel(axRaw, 'Speed (m/s)'); ylabel(axRaw, 'Range (m)');
    
    % Plot 2: MTI
    imagesc(axMTI, sGridMTI, rangeGridFiltered, 10*log10(magnitudeMTIMasked));
    axis(axMTI, 'xy'); xlabel(axMTI, 'Speed (m/s)');
    
    % Plot 3: CFAR (Showing only detections)
    imagesc(axCFAR, sGridMTI, rangeGridFiltered, detMapMasked);
    axis(axCFAR, 'xy'); xlabel(axCFAR, 'Speed (m/s)');
    
    % Update Alert UI (This part stays the same)
    set(hSquare, 'FaceColor', storedColor);
    set(hText, 'String', storedText);
    
    drawnow limitrate;
    
    % D. Target Extraction
    [row, col] = find(detMapMasked > 0);
    pts = [row, col]; 
    targets = [];
    
    if ~isempty(pts)
        try
            % Using min_points_cluster from the top
            idx = dbscan(pts, 5, min_points_cluster); 
            for k = 1:max(idx)
                clusterPts = pts(idx == k, :);
                centroidR = mean(clusterPts(:,1));
                targets = [targets; round(centroidR), 1, 1];
            end
        catch
            targets = target_extraction(detMapMasked, 5); 
        end
    end
    
    % E. State Machine (Diagnostic Mode)
    if ~isempty(targets)
        % Moving Target Detected
        lockedRangeIdx = targets(1, 1);
        realRange = rangeGridFiltered(lockedRangeIdx);
        
        if realRange <= 1.5,     currentZone = 3; % Red
        elseif realRange <= 3.5, currentZone = 2; % Yellow
        else,                    currentZone = 1; % Green
        end
        
        % For diagnostics, we'll bypass confirmation if it's 1
        lastDetectionTime = currentTime;
        isTargetLocked = true;
        
        if currentZone == 3
            storedColor = [1 0 0]; storedText = sprintf('DANGER: %.2fm', realRange);
        elseif currentZone == 2
            storedColor = [1 0.9 0]; storedText = sprintf('WARNING: %.2fm', realRange);
        else
            storedColor = [0 0.8 0]; storedText = sprintf('CLEAR: %.2fm', realRange);
        end
        
        % DEBUG PRINT
        fprintf('TARGET FOUND: %.2fm in Zone %d\n', realRange, currentZone);
        
    elseif isTargetLocked && (currentTime - lastDetectionTime) < 5.0
        % Stationary Path Logic
        [~, zeroDopplerIdx] = min(abs(sGridRaw)); 
        currentRawNoiseFloor = mean(magnitudeRaw(:));
        currentRawEnergy = magnitudeRaw(lockedRangeIdx, zeroDopplerIdx);
        
        if currentRawEnergy > (currentRawNoiseFloor * 1.2)
            if ~contains(storedText, '(STATIONARY)')
                storedText = [storedText, ' (STATIONARY)'];
            end
            lastDetectionTime = currentTime; 
        else
            isTargetLocked = false;
        end
    end
    
    % Global Reset if no signal for 2.5 seconds
    if (currentTime - lastDetectionTime) > 2.5
        storedColor = [0.8 0.8 0.8]; 
        storedText = 'Searching...';
        isTargetLocked = false;
    end
    
    % Update UI
    set(hSquare, 'FaceColor', storedColor);
    set(hText, 'String', storedText);
    drawnow limitrate; 
end

cleanupAntenna(rx, tx, bf, bf_TDD);