clear; close all; clc; 

dataFolder = ['C:\Users\leeju\Downloads\March New Data Set\Towards_Radar_10m'];
%Parameters
fc = 10e9; 
fs = 1e6; 
rampbandwidth = 2000e6; 
prf = 2000; 
tSweep = 5e-4; 
sweepSlope = rampbandwidth / tSweep;
applyMTI = true; 
calRange=2;
minRange = 0; 
maxRange = 10;
tCapture=30;
nPulses = 512;
mtiKernel = [1 -1]; 
ncoeff = length(mtiKernel);
dr = 3e8 / (2 * rampbandwidth);

rdRaw = phased.RangeDopplerResponse(DopplerOutput="Speed",...
    OperatingFrequency=fc,SampleRate=fs,RangeMethod="FFT",...
    SweepSlope=sweepSlope,PRFSource="Property",PRF=prf);

rdFiltered = phased.RangeDopplerResponse('DopplerOutput', 'Speed', 'OperatingFrequency', fc, ...
    'SampleRate', fs, 'RangeMethod', 'FFT', 'SweepSlope', sweepSlope, 'PRFSource', 'Property', 'PRF', prf);

scopeRaw = phased.RangeDopplerScope('IQDataInput', false, 'Name', '1. Raw Response', 'Position', [10 510 500 400]);
scopeMTI = phased.RangeDopplerScope('IQDataInput', false, 'Name', '2. MTI Filtered', 'Position', [520 510 500 400]);
scopeCFAR = phased.RangeDopplerScope('IQDataInput', false, 'Name', '3. CFAR Detections', 'Position', [520 50 500 400]);

fileList = dir(fullfile(dataFolder, 'data_frame_*.mat'));
numFiles = length(fileList);
frameIdx=1;
for i = 1:numFiles
    thisFile = fullfile(dataFolder, fileList(i).name);
    loadedStruct = load(thisFile);
    if ~isfield(loadedStruct, 'rawData'), continue; end
    rawData = loadedStruct.rawData;
    
    % 1. Raw Scope
    [respRaw, rGridRaw, sGridRaw] = rdRaw(rawData);
    rGridRaw=rGridRaw-calRange;
    maskRaw = rGridRaw >= minRange & rGridRaw <= maxRange;
    scopeRaw(abs(respRaw(maskRaw,:)), rGridRaw(maskRaw), sGridRaw);
   
    % 2. MTI Scope
    mtiData = filter(mtiKernel, 1, rawData, [], 2);
    mtiData = mtiData(:, ncoeff:end);
    [respMTI, rGridMTI, sGridMTI] = rdFiltered(mtiData);
    rGridMTI=rGridMTI-calRange;
    maskMTI = rGridMTI >= minRange & rGridMTI <= maxRange;
    magnitudeMTI = abs(respMTI(maskMTI, :));
    rangeMTI=rGridMTI(maskMTI);
    scopeMTI(magnitudeMTI, rangeMTI, sGridMTI);
    
    % 3. CFAR Scope
    [detMap, noiseMap] = cfar_filter(magnitudeMTI, 10, 2, 20, 10, 5e-3);
    
    detMapPlot = detMap;
    detMapPlot(detMapPlot <= 0) = 1;
    scopeCFAR(detMapPlot, rGridMTI(maskMTI), sGridMTI);

    targets = target_extraction(detMap, noiseMap, 12);
    if ~isempty(targets)
        fprintf('\n--- Frame %d: Found %d targets ---\n', frameIdx, size(targets,1));
        for k = 1:size(targets,1)
            r_idx = targets(k,1); % Range Index
            
            % Physical distance=Index * Resolution
            physical_range = r_idx*dr; 
            
            fprintf('   Target %d: Range Index = %d | Distance = %.2f meters\n', ...
                k, r_idx, physical_range);
        end
    end
    frameIdx=frameIdx+1;
    drawnow;
end