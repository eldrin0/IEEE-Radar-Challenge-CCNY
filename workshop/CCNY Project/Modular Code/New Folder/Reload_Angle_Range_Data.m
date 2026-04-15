clear; close all; clc;

% --- Playback Configuration ---
dataFolder = 'C:\Users\leeju\Documents\HeatMap Data'; % Ensure this matches your save location
playbackDelay = 0.05; % Delay in seconds between frames (adjust for playback speed)

% --- CFAR Filter Parameters ---
% Must match the parameters used during capture for accurate reproduction
guard_rng = 2;    
train_rng = 4;    
guard_az  = 1;    
train_az  = 2;    
Pfa       = 1e-3; 
% ------------------------------

% 1. Load Background Data (Need axes and mask)
bgFilename = fullfile(dataFolder, 'Background_Data.mat');
if ~isfile(bgFilename)
    error('Could not find Background_Data.mat in %s. Cannot proceed without axes data.', dataFolder);
end
disp('Loading Background Data...');
load(bgFilename, 'backgroundMask_dB', 'validRanges', 'azAngles');

% 2. Find all Raw Heatmap Data files
filePattern = fullfile(dataFolder, 'Raw_Heatmap_Data_*.mat');
dataFiles = dir(filePattern);
numFiles = length(dataFiles);

if numFiles == 0
    error('No Raw_Heatmap_Data files found in %s.', dataFolder);
end
fprintf('Found %d data frames. Starting playback...\n', numFiles);

% =========================================================
%                  SINGLE MAXIMIZED FIGURE
% =========================================================
% Initialize with a dummy matrix of zeros for the initial setup
dummyMatrix = zeros(length(validRanges), length(azAngles));

fMain = figure('Name', 'Radar Heatmap Playback', 'Color', 'w', 'WindowState', 'maximized');
t = tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

% Plot 1: Raw Heatmap
ax1 = nexttile(t);
hImage1 = imagesc(ax1, azAngles, validRanges, dummyMatrix);
set(ax1, 'YDir', 'normal'); 
colormap(ax1, 'jet');
cbar1 = colorbar(ax1);
cbar1.Label.String = 'Absolute Power (dB)';
clim(ax1, [30 90]);
xlabel(ax1, 'Azimuth Angle (Degrees)', 'FontWeight', 'bold');
ylabel(ax1, 'Range (Meters)', 'FontWeight', 'bold');
title(ax1, 'Raw 8-Element Heatmap', 'FontSize', 12);
peakText1 = text(ax1, 0.02, 0.95, 'Peak: -- °', 'Units', 'normalized', ...
    'FontSize', 12, 'FontWeight', 'bold', 'Color', 'white');

% Plot 2: Difference Heatmap
ax2 = nexttile(t);
hImage2 = imagesc(ax2, azAngles, validRanges, dummyMatrix);
set(ax2, 'YDir', 'normal'); 
colormap(ax2, 'jet');
cbar2 = colorbar(ax2);
cbar2.Label.String = 'Relative Power (dB)';
clim(ax2, 'auto'); 
xlabel(ax2, 'Azimuth Angle (Degrees)', 'FontWeight', 'bold');
ylabel(ax2, 'Range (Meters)', 'FontWeight', 'bold');
title(ax2, 'Difference Heatmap (Data - Mask)', 'FontSize', 12);
peakText2 = text(ax2, 0.02, 0.95, 'Peak: -- dB | R: -- m | Az: -- °', 'Units', 'normalized', ...
    'FontSize', 12, 'FontWeight', 'bold', 'Color', 'white');

% Plot 3: CFAR Filtered Heatmap
ax3 = nexttile(t);
hImage3 = imagesc(ax3, azAngles, validRanges, dummyMatrix);
set(ax3, 'YDir', 'normal'); 
colormap(ax3, 'jet');
cbar3 = colorbar(ax3);
cbar3.Label.String = 'Absolute Power (dB)';
clim(ax3, 'auto'); 
xlabel(ax3, 'Azimuth Angle (Degrees)', 'FontWeight', 'bold');
ylabel(ax3, 'Range (Meters)', 'FontWeight', 'bold');
title(ax3, 'CFAR Detection Map', 'FontSize', 12);
peakText3 = text(ax3, 0.02, 0.95, 'Peak: -- dB | R: -- m | Az: -- °', 'Units', 'normalized', ...
    'FontSize', 12, 'FontWeight', 'bold', 'Color', 'white');

% =========================================================
%                       PLAYBACK LOOP
% =========================================================

for k = 1:numFiles
    % Load current frame
    currentFile = fullfile(dataFiles(k).folder, dataFiles(k).name);
    load(currentFile, 'heatmapMatrix', 'heatmap_dB');
    
    % --- Calculate Difference ---
    difference_dB = heatmap_dB - backgroundMask_dB;
    difference_dB(difference_dB < 0) = 0;
    
    [maxDiffVal, linearIdx] = max(difference_dB(:));
    [maxRow, maxCol] = ind2sub(size(difference_dB), linearIdx);
    peakRange = validRanges(maxRow);
    peakAngleDiff = azAngles(maxCol);
    
    % --- Apply CFAR Filter ---
    [detection_map, ~] = cfar_filter(heatmapMatrix, guard_az, guard_rng, train_az, train_rng, Pfa);
    
    cfar_mask = detection_map > 0;
    currentMin = min(heatmap_dB(:));
    cfar_heatmap_dB = ones(size(heatmap_dB)) * currentMin; 
    cfar_heatmap_dB(cfar_mask) = heatmap_dB(cfar_mask);
    
    [maxCfarVal, cfarIdx] = max(cfar_heatmap_dB(:));
    if any(cfar_mask(:))
        [maxRowCfar, maxColCfar] = ind2sub(size(cfar_heatmap_dB), cfarIdx);
        peakRangeCfar = validRanges(maxRowCfar);
        peakAngleCfar = azAngles(maxColCfar);
        cfarText = sprintf('Peak: %5.1f dB | R: %4.2f m | Az: %3d°', maxCfarVal, peakRangeCfar, peakAngleCfar);
    else
        cfarText = 'No Targets Detected';
    end

    % --- Update Plots ---
    [~, maxColIdx1] = max(max(heatmap_dB)); 
    
    % Update Figure 1
    set(hImage1, 'CData', heatmap_dB);
    set(peakText1, 'String', sprintf('Peak: %d°', azAngles(maxColIdx1)));
    
    % Update Figure 2
    set(hImage2, 'CData', difference_dB);
    set(peakText2, 'String', sprintf('Peak: %5.1f dB | R: %4.2f m | Az: %3d°', maxDiffVal, peakRange, peakAngleDiff));
    
    % Update Figure 3
    set(hImage3, 'CData', cfar_heatmap_dB);
    set(peakText3, 'String', cfarText);
    
    % Update main title to show progress
    sgtitle(fMain, sprintf('Playback: Frame %d of %d', k, numFiles), 'FontSize', 14, 'FontWeight', 'bold');

    drawnow;
    
    % Add a short pause to control playback speed
    pause(playbackDelay);
    
    % Allow early exit if the user closes the figure
    if ~isvalid(fMain)
        disp('Playback terminated by user.');
        break;
    end
end

if isvalid(fMain)
    disp('Playback Complete.');
end

% =========================================================================
%                             FUNCTIONS
% =========================================================================

function [detection_map, noise_map] = cfar_filter(mag_data, guard_r, guard_d, train_r, train_d, Pfa)
% APPLY_CFAR  Applies CA-CFAR (Cell Averaging) detection on a Range-Doppler/Azimuth map.
    [Nd, Nr] = size(mag_data);
    
    win_d = guard_d + train_d;
    win_r = guard_r + train_r;
    
    % Total number of training cells
    N_train = (2*train_d + 2*guard_d + 1) * (2*train_r + 2*guard_r + 1) ...
            - (2*guard_d + 1)*(2*guard_r + 1);
            
    if N_train <= 0
        error('Number of training cells must be positive.');
    end
    
    % CA-CFAR threshold scaling factor
    alpha = N_train * (Pfa^(-1/N_train) - 1);
    kernel = ones(2*win_d + 1, 2*win_r + 1);
    
    center_d = win_d + 1;
    center_r = win_r + 1;
    kernel(center_d-guard_d : center_d+guard_d, center_r-guard_r : center_r+guard_r) = 0;
    
    kernel = kernel / N_train;
    noise_est_full = conv2(mag_data, kernel, 'same');
    noise_map = nan(Nd, Nr);
    
    valid_d = (1+win_d) : (Nd-win_d);
    valid_r = (1+win_r) : (Nr-win_r);
    
    noise_map(valid_d, valid_r) = noise_est_full(valid_d, valid_r);
    threshold_map = alpha * noise_map;
    
    cfar_mask = mag_data > threshold_map;
    detection_map = zeros(size(mag_data));
    detection_map(cfar_mask) = mag_data(cfar_mask);
end