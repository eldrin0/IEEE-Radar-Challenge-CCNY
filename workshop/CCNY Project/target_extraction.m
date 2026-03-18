function targets = target_extraction(detection_map, noise_map, max_targets)
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
        
        % --- THE FIX: Swap r_idx (Row/Range) and v_idx (Column/Doppler) ---
        [r_idx, v_idx] = ind2sub(size(temp_map), idx);

        % Look up the background noise level at this exact coordinate
        % Note: Matrices are always indexed (Row, Column), so it's (r_idx, v_idx)
        local_noise = noise_map(r_idx, v_idx);
        
        % Check if the peak is strong enough compared to the noise.
        if Valid_Target(peak_val, local_noise, 10) 
            % [Range Index, Doppler Index, Magnitude]
            targets = [targets; r_idx, v_idx, peak_val];
        end

        % Define a boundary of 10 pixels around the target we just found
        mask_w = 10;
        
        % Creates matrix around the valid target and obsucres it. 
        r1 = max(1, r_idx - mask_w);
        r2 = min(size(temp_map,1), r_idx + mask_w);
        v1 = max(1, v_idx - mask_w);
        v2 = min(size(temp_map,2), v_idx + mask_w);
        
        % Set everything inside that box to 0. 
        temp_map(r1:r2, v1:v2) = 0;
    end
end