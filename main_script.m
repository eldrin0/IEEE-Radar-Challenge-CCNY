clear; close all; clc;

% Define Radar Parameters
fc = 10e9; fs = 1e6; 
rampbandwidth = 2000e6; 
prf = 2000; tSweep = 5e-4; 
sweepSlope = rampbandwidth / tSweep;
calRange = 1.85; 
minRange = 0.2; maxRange = 6;
tCapture = 60;
nPulses = 512;
mtiKernel = [1 -1]; 
ncoeff = length(mtiKernel);

% Detection Logic & Stability Initialization
isTargetLocked = false;
lockedRangeIdx = 1;
storedColor = [0.8 0.8 0.8];
storedText = 'Searching...';
lastDetectionTime = -inf;
target_angle = 0;
azAngles = -60 : 5 : 60;

precomputedPhases = precompute_phases(fc, azAngles);

% Hardware & Scope Setup
[rx,tx,bf,bf_TDD] = setupLabRadar(fc,prf,nPulses,fs,rampbandwidth);

rdRaw = phased.RangeDopplerResponse(DopplerOutput="Speed",...
    OperatingFrequency=fc,SampleRate=fs,RangeMethod="FFT",...
    SweepSlope=sweepSlope,PRFSource="Property",PRF=prf);

rdFiltered = phased.RangeDopplerResponse('DopplerOutput', 'Speed', ...
    'OperatingFrequency', fc, 'SampleRate', fs, 'RangeMethod', 'FFT', ...
    'SweepSlope', sweepSlope, 'PRFSource', 'Property', 'PRF', prf);

% UI Setup
[hMainFig, axRaw, axMTI, axCFAR, axTarget, hSquare, hText] = setup_ui();
hImgRaw = []; hImgMTI = []; hImgCFAR = [];

% Main Processing Loop
amp = 0.9 * 2^15;
txWaveform = amp * ones(rx.SamplesPerFrame, 2);
tx(txWaveform);
t_start = tic;
steer(target_angle, bf, azAngles, precomputedPhases);

while toc(t_start) < tCapture
    % Capture Data
    data = capture_data(rx, bf, bf_TDD);
   
    % Raw Response & MTI Filtering
    [magnitudeRaw, sGridRaw, rangeGridFiltered, detMapMasked, ...
        magnitudeMTIMasked, sGridMTI] = process_signals(data, rdRaw, ...
        rdFiltered, calRange, minRange, maxRange, mtiKernel, ncoeff);
    
    % Plot Raw, MTI & CFAR Responses
    [hImgRaw, hImgMTI, hImgCFAR] = update_dashboard(hImgRaw, hImgMTI, ...
        hImgCFAR, axRaw, axMTI, axCFAR, sGridRaw, rangeGridFiltered, ...
        magnitudeRaw, sGridMTI, magnitudeMTIMasked, detMapMasked, ...
        hSquare, storedColor, hText, storedText);
    
    % Extract Targets
    targets = extract_targets(detMapMasked);
    
    % Updated detection state (Red, Yellow, Green)
    [lockedRangeIdx, lastDetectionTime, isTargetLocked, storedColor, ...
        storedText] = update_state(targets, lockedRangeIdx, ...
                      rangeGridFiltered, currentTime, lastDetectionTime, ...
                      isTargetLocked, sGridRaw, magnitudeRaw, ...
                      storedColor, storedText);
end

cleanupAntenna(rx, tx, bf, bf_TDD);