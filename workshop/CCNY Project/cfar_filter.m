function [detection_map, noise_map] = cfar_filter(mag_data, guard_r, guard_d, train_r, train_d, Pfa)
% APPLY_CFAR  Applies CA-CFAR (Cell Averaging) detection on a Range-Doppler map.

    [Nd, Nr] = size(mag_data);
    
    win_d = guard_d + train_d;
    win_r = guard_r + train_r;

    % Total number of training cells
    N_train = (2*train_d + 2*guard_d + 1) * (2*train_r + 2*guard_r + 1) ...
            - (2*guard_d + 1)*(2*guard_r + 1);

    if N_train <= 0
        error('Number of training cells must be positive.');
    end

    % CA-CFAR threshold scaling factor
    alpha = N_train * (Pfa^(-1/N_train) - 1);

    kernel = ones(2*win_d + 1, 2*win_r + 1);
    
    center_d = win_d + 1;
    center_r = win_r + 1;
    kernel(center_d-guard_d : center_d+guard_d, center_r-guard_r : center_r+guard_r) = 0;
    
    kernel = kernel / N_train;

    noise_est_full = conv2(mag_data, kernel, 'same');

    noise_map = nan(Nd, Nr);
    valid_d = (1+win_d) : (Nd-win_d);
    valid_r = (1+win_r) : (Nr-win_r);
    
    noise_map(valid_d, valid_r) = noise_est_full(valid_d, valid_r);

    threshold_map = alpha * noise_map;
    cfar_mask = mag_data > threshold_map;

    detection_map = zeros(size(mag_data));
    detection_map(cfar_mask) = mag_data(cfar_mask);
end


