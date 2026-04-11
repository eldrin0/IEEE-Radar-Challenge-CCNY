function fns = radar_setup()
    % Return handles to the processing functions below
    fns.initialize = @initialize_radar_system;
    fns.collect    = @collect_radar_frame;
    fns.processRaw = @process_raw_response;
    fns.cleanup    = @cleanup_radar_system;
end

function [ui, hw, dsp, grid] = initialize_radar_system(p)
    % 1. Setup UI via the dashboard library
    gui_lib = radar_dashboard(); 
    ui = gui_lib.setup();
    
    % 2. Setup Hardware
    [hw.rx, hw.tx, hw.bf, hw.bf_TDD] = setupLabRadar(p.fc, p.prf, p.nPulses, p.fs, p.bw);
    
    % 3. Setup DSP Objects (FIX: Independent initialization)
    slope = p.bw / p.tSweep;
    
    dsp.rdRaw = phased.RangeDopplerResponse('DopplerOutput','Speed',...
        'OperatingFrequency',p.fc,'SampleRate',p.fs,'RangeMethod','FFT',...
        'SweepSlope',slope,'PRFSource','Property','PRF',p.prf);
        
    dsp.rdFiltered = phased.RangeDopplerResponse('DopplerOutput','Speed',...
        'OperatingFrequency',p.fc,'SampleRate',p.fs,'RangeMethod','FFT',...
        'SweepSlope',slope,'PRFSource','Property','PRF',p.prf);
    
    % 4. Pre-Calculations & Masking (FIX for checkPRF error)
    maxSamples = floor(p.fs / p.prf);
    [~, rGridTemp, sGrid] = dsp.rdRaw(zeros(maxSamples, p.nPulses));
    
    rGridCal = rGridTemp - p.calRange;
    grid.mask = rGridCal >= p.minRange & rGridCal <= p.maxRange;
    grid.sGrid = sGrid;
    grid.rGridFull = rGridCal;
    
    % 5. Start Transmitter
    txWaveform = (0.9 * 2^15) * ones(hw.rx.SamplesPerFrame, 2);
    hw.tx(txWaveform);
end

function data = collect_radar_frame(hw, p)
    try
        hw.bf.Burst = false; hw.bf.Burst = true; hw.bf.Burst = false;
        raw = hw.rx();
        dataFull = arrangePulseData(raw, hw.rx, hw.bf, hw.bf_TDD);
        
        % Truncate to ensure DSP compatibility
        maxS = floor(p.fs / p.prf);
        data = dataFull(1:maxS, :);
    catch
        data = [];
    end
end

function [magRaw, rGrid] = process_raw_response(data, rdObj, grid)
    [respRaw, ~, ~] = rdObj(data);
    magRaw = abs(respRaw(grid.mask, :));
    rGrid = grid.rGridFull(grid.mask);
end

function cleanup_radar_system(hw)
    cleanupAntenna(hw.rx, hw.tx, hw.bf, hw.bf_TDD);
    disp('Radar Session Ended Safely.');
end