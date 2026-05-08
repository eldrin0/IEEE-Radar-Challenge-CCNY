function targets = extract_targets(detMapMasked)
    [row, col] = find(detMapMasked > 0);
    pts = [row, col]; targets = [];
    if ~isempty(pts)
        try
            idx = dbscan(pts, 5, 1); 
            for k = 1:max(idx)
                clusterPts = pts(idx == k, :);
                centroidR = mean(clusterPts(:,1));
                targets = [targets; round(centroidR), 1, 1];
            end
        catch
            temp_map = detMapMasked; 
            targets = [];
            for t = 1:5
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
    end
end