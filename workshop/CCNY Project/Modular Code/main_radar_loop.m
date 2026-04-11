clear all; close all; clc;



% --- 1. Load Libraries ---
gui = radar_dashboard();  
eng = radar_setup(); 



% --- 2. Configuration & Initialization ---
params = struct('fc', 10e9, 'fs', 1e6, 'prf', 2000, 'bw', 2000e6, 'tSweep', 5e-4, 'nPulses', 512, 'calRange', 1.85, 'minRange', 0.2, 'maxRange', 10);

[ui, hw, dsp, grid] = eng.initialize(params);




% --- 3. Main Loop ---
t_start = tic;
tCapture = 600;
clear detection_state; 

while toc(t_start) < tCapture
    currentTime = toc(t_start);
    
    % Get Data
    data = eng.collect(hw, params);
    if isempty(data), continue; end
    
    % Signal Processing
    [magRaw, rGrid] = eng.processRaw(data, dsp.rdRaw, grid);
    [targets, detMap, magMTI] = process_radar_data(data, dsp.rdFiltered, [1 -1], grid.mask, 5.0e-5);
    
    % Decision Logic for Safety State
    [color, text] = detection_state(targets, magRaw, rGrid, currentTime);
    
    % Update Dashboard
    gui.update(ui, grid.sGrid, rGrid, magRaw, magMTI, detMap, color, text);
    
    drawnow limitrate;
end




% --- 4. Shutdown ---
eng.cleanup(hw);