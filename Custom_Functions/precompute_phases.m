function precomputedPhases = precompute_phases(fc, azAngles)
    numElements = 8;
    elementSpacing = 0.014;
    array = phased.ULA('NumElements', numElements, 'ElementSpacing', elementSpacing);
    elementLocation = array.getElementPosition();
    load('CalibrationWeights.mat', 'calibrationweights');
    c = physconst('LightSpeed');
    lambda = c / fc;
    
    disp('Pre-computing phase lookup tables...');
    precomputedPhases = zeros(length(azAngles), 8);
    for i = 1:length(azAngles)
        steerweights = steervec(elementLocation / lambda, [azAngles(i); 0]);
        steer_cols = [steerweights(1:4) steerweights(5:8)];
        analogWeights = analogWeightsCalAdjustment(steer_cols, calibrationweights.AnalogWeights);
        analog_flat = [analogWeights(:,1); analogWeights(:,2)].';
        precomputedPhases(i, :) = mod(rad2deg(angle(analog_flat)), 360);
    end
end