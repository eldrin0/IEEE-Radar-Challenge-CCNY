function steer(target_angle, bf, azAngles, precomputedPhases)
    idx = find(azAngles == target_angle, 1);
    phase_deg = precomputedPhases(idx, :);
    bf.RxPhase = phase_deg;
    bf.LatchRxSettings();
end