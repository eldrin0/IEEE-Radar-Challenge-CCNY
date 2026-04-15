%% Search & Track Radar System
clear; close all; clc;
warning('off','MATLAB:system:ObsoleteSystemObjectMixin');

fc = 10e9;
prf = 2000;
fs = 1e6;
rampbandwidth = 2500e6;
c = physconst('LightSpeed');
lambda = c / fc;
elementSpacing = 0.014;
numElements = 8;
[rx, tx, bf, bf_TDD] = setupLabRadar(fc, prf, 1, fs, rampbandwidth);

tCapture = 120;
minRange = 0.2;
maxRange = 6;
calRange = 1.85;
tSweep = double(bf.FrequencyDeviationTime)/1e6;
sweepSlope = rampbandwidth/tSweep;

array = phased.ULA('NumElements', numElements, 'ElementSpacing', elementSpacing);
elementLocation = array.getElementPosition();
rr = phased.RangeResponse('RangeMethod', 'FFT', 'SweepSlope', sweepSlope, 'SampleRate', fs);

load('CalibrationWeights.mat', 'calibrationweights');
elementEnable = [1 1; 1 1; 1 1; 1 1];
digitalsteer = digitalWeightsCalAdjustment([1;1], calibrationweights.DigitalWeights);

beginningAngle = -60;
endAngle = 60;
incrementalAngle = 5;

azAngles = beginningAngle : incrementalAngle : endAngle; 
nAngles = length(azAngles);

disp('Pre-computing phase lookup tables...');
precomputedPhases = zeros(nAngles, 8);
for i = 1:nAngles
    steerAngle = azAngles(i);
    steerweights = steervec(elementLocation / lambda, [steerAngle; 0]);
    steer_cols = [steerweights(1:4) steerweights(5:8)];
    analogsteer = analogWeightsCalAdjustment(steer_cols, calibrationweights.AnalogWeights);
    analogWeights = analogsteer .* elementEnable;
    
    analog_flat = [analogWeights(:,1); analogWeights(:,2)].';
    precomputedPhases(i, :) = mod(rad2deg(angle(analog_flat)), 360);
end

amp = 0.9 * 2^15;
txWaveform = amp*ones(rx.SamplesPerFrame,2);
tx(txWaveform); 
bf.Burst = true; 
testData = rx();
bf.Burst = false;

testDataArranged = arrangePulseData(testData, rx, bf, bf_TDD);
[~, range] = rr(testDataArranged(:,1));
range = range - calRange;
keepRange = range >= minRange & range <= maxRange;
validRanges = range(keepRange); 

heatmapMatrix = zeros(length(validRanges), nAngles);

f1 = figure('Name', 'Live Calibrated Heatmap', 'Color', 'w', 'Position', [100, 100, 700, 600]);
ax1 = axes(f1);
hImage1 = imagesc(ax1, azAngles, validRanges, heatmapMatrix);
set(ax1, 'YDir', 'normal'); % Keep range 0 at the bottom
colormap(ax1, 'jet');
cbar1 = colorbar(ax1);
cbar1.Label.String = 'Absolute Power (dB)';
clim(ax1, [35 100]); 
xlabel(ax1, 'Azimuth Angle (Degrees)', 'FontWeight', 'bold');
ylabel(ax1, 'Range (Meters)', 'FontWeight', 'bold');
title(ax1, 'Raw 8-Element Heatmap', 'FontSize', 12);
peakText1 = text(ax1, 0.02, 0.95, 'Peak: -- °', 'Units', 'normalized', ...
    'FontSize', 12, 'FontWeight', 'bold', 'Color', 'white');

f2 = figure('Name', 'Difference Heatmap (Target Detection)', 'Color', 'w', 'Position', [820, 100, 700, 600]);
ax2 = axes(f2);
hImage2 = imagesc(ax2, azAngles, validRanges, heatmapMatrix);
set(ax2, 'YDir', 'normal'); 
colormap(ax2, 'jet');
cbar2 = colorbar(ax2);
cbar2.Label.String = 'Relative Power (dB)';
clim(ax2, [0 25]); % Scaled for difference
xlabel(ax2, 'Azimuth Angle (Degrees)', 'FontWeight', 'bold');
ylabel(ax2, 'Range (Meters)', 'FontWeight', 'bold');
title(ax2, 'Difference Heatmap (Live - Mask)', 'FontSize', 12);
peakText2 = text(ax2, 0.02, 0.95, 'Peak: -- dB | R: -- m | Az: -- °', 'Units', 'normalized', ...
    'FontSize', 12, 'FontWeight', 'bold', 'Color', 'white');

captureTransmitWaveform(rx, tx, bf);
disp('---------------------------------------------------------------------------------');
disp(' Total Sweep | Analog Setup | Hardware Capture | Digital Proc. |  Plot Update  ');
disp('---------------------------------------------------------------------------------');

tTotal = tic; 
isFirstSweep = true; 
backgroundMask_dB = zeros(length(validRanges), nAngles);

while toc(tTotal) < tCapture
    sweepStart = tic;

    [heatmap_dB, heatmapMatrix, time_analog, time_capture, time_digital] = getHeatmapData(nAngles, precomputedPhases, validRanges, rx, bf, bf_TDD, rr, keepRange, digitalsteer);
    
    tP = tic;
    
    if isFirstSweep
        backgroundMask_dB = heatmap_dB;
        disp('Initial background mask captured!');
        isFirstSweep = false;
    end
    
    difference_dB = heatmap_dB - backgroundMask_dB;
    difference_dB(difference_dB < 0) = 0;
    
    [maxDiffVal, linearIdx] = max(difference_dB(:));
    [maxRow, maxCol] = ind2sub(size(difference_dB), linearIdx);
    peakRange = validRanges(maxRow);
    peakAngleDiff = azAngles(maxCol);
    
    [~, maxColIdx1] = max(max(heatmap_dB)); 
    set(hImage1, 'CData', heatmap_dB);
    set(peakText1, 'String', sprintf('Peak: %d°', azAngles(maxColIdx1)));
    
    set(hImage2, 'CData', difference_dB);
    set(peakText2, 'String', sprintf('Peak: %5.1f dB | R: %4.2f m | Az: %3d°', maxDiffVal, peakRange, peakAngleDiff));
    
    drawnow;
    time_plot = toc(tP);
    
    time_total_sweep = toc(sweepStart);
    
    fprintf('   %5.3f s   |   %5.3f s    |     %5.3f s      |   %5.3f s    |   %5.3f s\n', ...
        time_total_sweep, time_analog, time_capture, time_digital, time_plot);
end

disp('Capture Complete.');
cleanupAntenna(rx, tx, bf, bf_TDD);

function [heatmap_dB, heatmapMatrix, tA_total, tC_total, tD_total] = getHeatmapData(nAngles, precomputedPhases, validRanges, rx, bf, bf_TDD, rr, keepRange, digitalsteer)
    
    heatmapMatrix = zeros(length(validRanges), nAngles);
    
    tA_total = 0; 
    tC_total = 0; 
    tD_total = 0;
    
    for i = 1:nAngles
        tA = tic;
        bf.RxPhase = precomputedPhases(i, :);
        bf.LatchRxSettings();
        tA_total = tA_total + toc(tA);
        
        tC = tic;
        bf.Burst = true; 
        data = rx();
        bf.Burst = false;
        tC_total = tC_total + toc(tC);
        
        tD = tic;
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
        
        heatmapMatrix(:, i) = abs(resp(keepRange));
        tD_total = tD_total + toc(tD);
    end
    
    heatmap_dB = 20 * log10(heatmapMatrix + 1e-9); 
end