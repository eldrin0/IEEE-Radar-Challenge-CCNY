function targets = target_extraction(detection_map, max_targets)
    temp_map = detection_map; 
    targets = [];
    for t = 1:max_targets
        [peak_val, idx] = max(temp_map(:));
        if peak_val == 0, break; end
        [r_idx, v_idx] = ind2sub(size(temp_map), idx);
        targets = [targets; r_idx, v_idx, peak_val];
        mask_w = 8;
        r1 = max(1, r_idx-mask_w); r2 = min(size(temp_map,1), r_idx+mask_w);
        v1 = max(1, v_idx-mask_w); v2 = min(size(temp_map,2), v_idx+mask_w);
        temp_map(r1:r2, v1:v2) = 0;
    end
end