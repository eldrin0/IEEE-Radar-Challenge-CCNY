function steer(target_angle, hw, azAngles, precomputedPhases)
    % Steers the radar by finding the exact precomputed phase values
    
    idx = find(azAngles == target_angle, 1);
    phase_deg = precomputedPhases(idx, :);
    hw.bf.RxPhase = phase_deg;
    hw.bf.LatchRxSettings();
    
    fprintf('Radar beam steered to %d°\n', target_angle);
end