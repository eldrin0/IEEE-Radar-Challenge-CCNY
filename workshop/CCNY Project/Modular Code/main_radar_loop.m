clear all; close all; clc;


% --- 1. Load Libraries ---
gui = radar_dashboard();  
eng = radar_setup(); 


% --- 2. Configuration & Initialization ---
params = struct('fc', 10e9, 'fs', 1e6, 'prf', 2000, 'bw', 2000e6, 'tSweep', 5e-4, 'nPulses', 512, 'calRange', 1.85, 'minRange', 0.2, 'maxRange', 10);

tic
[ui, hw, dsp, grid] = eng.initialize(params);
endTime = toc;
fprintf('Initialization Time: %.4f seconds\n', endTime);



% --- 3. Main Loop ---
t_start = tic;
tCapture = 600;
clear detection_state; 

fprintf('Data Collection  |  Signal Processing  |  Dashboard Update  |  Total\n');

while toc(t_start) < tCapture
    currentTime = toc(t_start);
    
    % Get Data
    tic
    data = eng.collect(hw, params);
    if isempty(data), continue; end
    endTime1 = toc;

    % Signal Processing
    tic
    [magRaw, rGrid] = eng.processRaw(data, dsp.rdRaw, grid);
    endTime2 = toc;


    [targets, detMap, magMTI] = process_radar_data(data, dsp.rdFiltered, [1 -1], grid.mask, 5.0e-5);
    
    % Decision Logic for Safety State
    [color, text] = detection_state(targets, magRaw, rGrid, currentTime);
    
    % Update Dashboard
    tic
    gui.update(ui, grid.sGrid, rGrid, magRaw, magMTI, detMap, color, text);
    endTime3 = toc;
    totalTime = endTime1 + endTime2 + endTime3;
    % Log processing times for performance analysis
    fprintf('D_Coll%.4f               %.4f              %.4f                %.4f', endTime1, endTime2, endTime3, totalTime);
    fprintf('\n');
    drawnow limitrate;
end




% --- 4. Shutdown ---
eng.cleanup(hw);