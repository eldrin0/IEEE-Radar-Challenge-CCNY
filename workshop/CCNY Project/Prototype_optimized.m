cleanupAntenna(rx,tx,bf,bf_TDD);

clear;
close all;
warning('off','MATLAB:system:ObsoleteSystemObjectMixin');


load('CalibrationWeights.mat','calibrationweights');
fc = 10e9;
prf = 4000;
nPulses = 512;
calrange = 2;
fs = 1e6;
rampbandwidth = 500e6;
tCapture = 30;
displayRate = 4; % Update the plot every 4th frame


[rx,tx,bf,bf_TDD] = setupLabRadar(fc, prf, nPulses, fs, rampbandwidth);

mti = [1 -2 1];
ncoeff = length(mti);
applyFilter = false;
minIntensity = 30;
maxIntensity = 100;

minRange = 0;
maxRange = 10;


tSweep = double(bf.FrequencyDeviationTime)/1e6;
sweepSlope = rampbandwidth/tSweep;


rd = phased.RangeDopplerResponse(DopplerOutput="Speed",...
    OperatingFrequency=fc,SampleRate=fs,RangeMethod="FFT",...
    SweepSlope=sweepSlope,PRFSource="Property",PRF=prf);


rr = phased.RangeResponse(RangeMethod="FFT",SweepSlope=sweepSlope,SampleRate=fs);


scope = phased.RangeDopplerScope(IQDataInput=false);
scope.Position = [50 50 800 450];

calRange=0;

% Create axes to plot range response and beat frequency on the same figure
f = figure('Name', 'Real Time Monitoring');
rangeAx = axes(f); hold(rangeAx,"on");
freqAx = axes(f); hold(freqAx,"on");
rangeline = plot(rangeAx, 0, 0, 'LineWidth', 1.5);
freqline = plot(freqAx, 0, 0, 'Color', 'r', 'LineStyle', '--');
rangeAx.Position = freqAx.Position;
xlabel(rangeAx,'Range (m)'); ylabel(rangeAx,'Magnitude');
freqAx.XAxisLocation = "top"; freqAx.YTick = []; xlabel(freqAx,'Beat Freq (KHz)');


maxFrames = ceil(tCapture * 20);
dataBuffer = complex(zeros(nPulses, 512, maxFrames, 'single'));

captureTransmitWaveform(rx,tx,bf);

tStart = tic;
frameIdx = 1;

%Name the output folder
%outputFolder = 'C:\Users\eldri\Downloads\Phaser-Control-with-MATLAB-main\Phaser-Control-with-MATLAB-main\Collected_Data\Chamber\EmptyRoom';
%if ~exist(outputFolder, 'dir')
%    mkdir(outputFolder);
%end

while toc(tStart) < tCapture
    rawFrame = captureTransmitWaveform(rx, tx, bf);
    data = arrangePulseData(rawFrame, rx, bf, bf_TDD);
    
    dataBuffer(:,:,frameIdx) = single(data);

    %fname = sprintf('data_frame_%04d.mat', frameIdx);
    %fullPath = fullfile(outputFolder, fname);
    displayData = data;
    if applyFilter
            displayData = filter(mti, 1, data, [], 2);
            displayData = displayData(:, ncoeff:end);
    end
    
    if mod(frameIdx, displayRate) == 0
    %rawData = data;
    %save(fullPath, 'rawData');
        [resp, range, speed] = rd(displayData);
        range = range - calrange;
        keepRange = range >= minRange & range <= maxRange;
        keepRangeData = resp(keepRange,:);
        scope(keepRangeData,range(keepRange),speed);
    

        % Get the range-Response for just one pulse
        [resp,range] = rr(data(:,1));
        % Generate the beat frequencies corresponding to ranges
        bfreq = linspace(-rr.SampleRate/2,rr.SampleRate/2,length(range));
        % Adjust the range for calibration range
        range = range - calRange;
        % Get range of interest
        keepRange = range >= minRange & range <= maxRange;

        % Plot the response on the range and freq axes
        set(rangeline, 'XData', range(keepRange), 'YData', abs(resp(keepRange)));
        set(freqline, 'XData', bfreq(keepRange)/1e3, 'YData', abs(resp(keepRange)));
        drawnow limitrate;
    end
    frameIdx = frameIdx + 1;
end


cleanupAntenna(rx,tx,bf,bf_TDD);