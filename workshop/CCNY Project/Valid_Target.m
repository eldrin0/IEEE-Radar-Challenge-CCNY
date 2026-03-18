function is_valid = Valid_Target(peak_value, noise_floor, min_snr_db)

    % --- Default ---
    if nargin < 3
        min_snr_db = 12.0;
    end

    if ~isscalar(peak_value) || peak_value < 0
        error('peak_value must be a non-negative scalar.');
    end
    if ~isscalar(noise_floor) || noise_floor < 0
        error('noise_floor must be a non-negative scalar.');
    end
    if ~isscalar(min_snr_db) || min_snr_db < 0
        error('min_snr_db must be a non-negative scalar.');
    end

    % --- Calculate SNR ---
    snr_linear = peak_value / (noise_floor + 1e-12);  % evita divisione per zero
    snr_db = 10 * log10(snr_linear);

    % --- Verify Output
    is_valid = snr_db >= min_snr_db;
end