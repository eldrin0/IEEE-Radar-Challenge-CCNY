clear; close all; clc;

%Parameters
fc = 10e9; 
fs = 1e6; 
rampbandwidth = 2000e6; 
prf = 2000; 
tSweep = 1e-3; 
sweepSlope = rampbandwidth / tSweep;
applyMTI = true; 
calRange=4;
minRange = 0; 
maxRange = 10;
tCapture=30;
nPulses = 512;
mtiKernel = [1 -2 1]; 
ncoeff = length(mtiKernel);
dr = 3e8 / (2 * rampbandwidth);

[rx,tx,bf,bf_TDD] = setupLabRadar(fc,prf,nPulses,fs,rampbandwidth);

Raw = captureTransmitWaveform(rx,tx,bf);
Data = arrangePulseData(Raw,rx,bf,bf_TDD);
[nSamples, nPulses] = size(Data);

rdRaw = phased.RangeDopplerResponse(DopplerOutput="Speed",...
    OperatingFrequency=fc,SampleRate=fs,RangeMethod="FFT",...
    SweepSlope=sweepSlope,PRFSource="Property",PRF=prf);

rdFiltered = phased.RangeDopplerResponse('DopplerOutput', 'Speed', 'OperatingFrequency', fc, ...
    'SampleRate', fs, 'RangeMethod', 'FFT', 'SweepSlope', sweepSlope, 'PRFSource', 'Property', 'PRF', prf);

scopeRaw = phased.RangeDopplerScope('IQDataInput', false, 'Name', '1. Raw Response', 'Position', [10 510 500 400]);
%scopeMTI = phased.RangeDopplerScope('IQDataInput', false, 'Name', '2. MTI Filtered', 'Position', [520 510 500 400]);
%scopeCFAR = phased.RangeDopplerScope('IQDataInput', false, 'Name', '3. CFAR Detections', 'Position', [520 50 500 400]);

captureTransmitWaveform(rx,tx,bf);

outputFolder = 'C:\Users\eldri\Downloads\NewDataSet\Moving_Away_Radar_4m';
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end
frameIdx = 1;
t = tic;
while toc(t) < tCapture
    data = captureTransmitWaveform(rx,tx,bf);
    data = arrangePulseData(data,rx,bf,bf_TDD);
    
    fname = sprintf('data_frame_%04d.mat', frameIdx);
    fullPath = fullfile(outputFolder, fname);
    
    rawData = data;
    save(fullPath, 'rawData');
    
    frameIdx = frameIdx + 1;
    % 1. Raw Scope
    [respRaw, rGridRaw, sGridRaw] = rdRaw(rawData);
    rGridRaw=rGridRaw-calRange;
    maskRaw = rGridRaw >= minRange & rGridRaw <= maxRange;
    scopeRaw(abs(respRaw(maskRaw,:)), rGridRaw(maskRaw), sGridRaw);
   
    % 2. MTI Scope
    mtiData = filter(mtiKernel, 1, rawData, [], 2);
    mtiData = mtiData(:, ncoeff:end);
    [respMTI, rGridMTI, sGridMTI] = rdFiltered(mtiData);
    maskMTI = rGridMTI >= minRange & rGridMTI <= maxRange;
    magnitudeMTI = abs(respMTI(maskMTI, :));
    %scopeMTI(magnitudeMTI, rGridMTI(maskMTI), sGridMTI);
    
    
    
    % 3. CFAR Scope
    [detMap, noiseMap] = apply_cfar(magnitudeMTI, 2, 2, 10, 10, 1e-4);
    %scopeCFAR(detMap, rGridMTI(maskMTI), sGridMTI);
    
    
    targets = extract_targets(detMap, noiseMap, 10);
    if ~isempty(targets)
        fprintf('\n--- Frame %d: Found %d targets ---\n', i, size(targets,1));
        for k = 1:size(targets,1)
            r_idx = targets(k,1); % Range Index
            
            % Physical distance=Index * Resolution
            physical_range = (r_idx - 1) * dr; 
            
            fprintf('   Target %d: Range Index = %d | Distance = %.2f meters\n', ...
                k, r_idx, physical_range);
        end
    end

    drawnow;
    pause(0.1); 
end
cleanupAntenna(rx,tx,bf,bf_TDD);

function targets = extract_targets(detection_map, noise_map, max_targets)
% EXTRACT_TARGETS  Extracts the strongest targets from a CFAR detection map.
%   Outputs a matrix where each row is a target: [Range Index, Doppler Index, Peak Value]

    % (setting pixels to zero) without destroying the original data.
    temp_map = detection_map; 
    
    % Initialize an empty array to store our valid targets
    targets = [];

    for t = 1:max_targets
        
        % Find the absolute highest value in the entire map, and get its linear index
        [peak_val, idx] = max(temp_map(:));
        
        % If the highest value is 0, the map is empty. Stop looking for targets.
        if peak_val == 0
            break; 
        end
        
        % Convert the 1D linear index back into 2D row (v_idx) and column (r_idx) coordinates
        [v_idx, r_idx] = ind2sub(size(temp_map), idx);

        % Look up the background noise level at this exact coordinate
        local_noise = noise_map(v_idx, r_idx);
        
        % Check if the peak is strong enough compared to the noise.
        if is_valid_detection(peak_val, local_noise, 10) 
            % If it passes, add a new row to our targets list:
            % [Column/Range Index, Row/Doppler Index, Magnitude]
            targets = [targets; r_idx, v_idx, peak_val];
        end

        % Define a boundary of 10 pixels around the target we just found
        mask_w = 10;
        
        % Calculate the top, bottom, left, and right edges of this boundary box.
        % The max() and min() functions ensure the box doesn't accidentally 
        % try to go outside the edges of the matrix.
        v1 = max(1, v_idx - mask_w);                  % Top edge
        v2 = min(size(temp_map,1), v_idx + mask_w);   % Bottom edge
        r1 = max(1, r_idx - mask_w);                  % Left edge
        r2 = min(size(temp_map,2), r_idx + mask_w);   % Right edge
        
        % Set everything inside that box to 0. 
        % This prevents the loop from finding the edges of this same target on the next run.
        temp_map(v1:v2, r1:r2) = 0;
    end
end

function [detection_map, noise_map] = apply_cfar(mag_data, guard_r, guard_d, train_r, train_d, Pfa)
% APPLY_CFAR  Applies CA-CFAR (Cell Averaging) detection on a Range-Doppler map.
%
%   [detection_map, noise_map] = APPLY_CFAR(mag_data, guard_r, guard_d, train_r, train_d, Pfa)
%
%   Inputs:
%       mag_data     - Magnitude map (e.g., |RD|), size [Nd x Nr]
%       guard_r      - Number of guard cells in range
%       guard_d      - Number of guard cells in Doppler
%       train_r      - Number of training cells in range
%       train_d      - Number of training cells in Doppler
%       Pfa          - Desired false alarm rate (e.g., 1e-5)
%
%   Outputs:
%       detection_map - Map with detections (same size as mag_data)
%       noise_map     - Estimated noise level at each cell

    if nargin < 6
        error('apply_cfar:MissingInputs', 'All six input arguments are required.');
    end
    if ~ismatrix(mag_data)
        error('mag_data must be a 2D matrix.');
    end
    if any([guard_r, guard_d, train_r, train_d] < 0) || ...
       ~all(fix([guard_r, guard_d, train_r, train_d]) == [guard_r, guard_d, train_r, train_d])
        error('Guard and training cell values must be non-negative integers.');
    end
    if ~isscalar(Pfa) || Pfa <= 0 || Pfa >= 1
        error('Pfa must be a scalar in the range (0, 1).');
    end

    [Nd, Nr] = size(mag_data);
    cfar_mask = false(Nd, Nr);
    noise_map = nan(Nd, Nr);  % Initialize noise threshold map

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

    % Optional masking of close range artifacts
    % mag_data(:, 1:25) = 0;

    % --- CFAR processing loop ---
    for i = 1+win_d : Nd-win_d
        for j = 1+win_r : Nr-win_r
            CUT = mag_data(i, j);

            % Extract local window
            local = mag_data(i-win_d:i+win_d, j-win_r:j+win_r);

            % Exclude guard cells and CUT
            local((win_d-guard_d+1):(win_d+guard_d+1), ...
                  (win_r-guard_r+1):(win_r+guard_r+1)) = NaN;

            noise_cells = local(~isnan(local));

            if isempty(noise_cells)
                continue;  % skip if not enough training cells
            end

            % Noise estimation (CA-CFAR)
            noise_est = mean(noise_cells);
            threshold = alpha * noise_est;

            noise_map(i, j) = noise_est;

            if CUT > threshold
                cfar_mask(i, j) = true;
            end
        end
    end

    % Output detection map
    detection_map = zeros(size(mag_data));
    detection_map(cfar_mask) = mag_data(cfar_mask);
end

function is_valid = is_valid_detection(peak_value, noise_floor, min_snr_db)
% IS_VALID_DETECTION  Valuta se una detection ha un SNR sufficiente
%
%   is_valid = IS_VALID_DETECTION(peak_value, noise_floor, min_snr_db)
%
%   Inputs:
%       peak_value   - Valore del picco rilevato
%       noise_floor  - Stima del rumore locale
%       min_snr_db   - Soglia minima di SNR in dB (default: 12 dB)
%
%   Output:
%       is_valid     - true se il picco ha SNR sufficiente, false altrimenti

    % --- Default ---
    if nargin < 3
        min_snr_db = 12.0;
    end

    % --- Sanificazione input ---
    if ~isscalar(peak_value) || peak_value < 0
        error('peak_value must be a non-negative scalar.');
    end
    if ~isscalar(noise_floor) || noise_floor < 0
        error('noise_floor must be a non-negative scalar.');
    end
    if ~isscalar(min_snr_db) || min_snr_db < 0
        error('min_snr_db must be a non-negative scalar.');
    end

    % --- Calcolo SNR ---
    snr_linear = peak_value / (noise_floor + 1e-12);  % evita divisione per zero
    snr_db = 10 * log10(snr_linear);

    % --- Verifica soglia ---
    is_valid = snr_db >= min_snr_db;
end
