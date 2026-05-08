function [magnitudeRaw, sGridRaw, rangeGridFiltered, detMapMasked, magnitudeMTIMasked, sGridMTI] = process_signals(data, rdRaw, rdFiltered, calRange, minRange, maxRange, mtiKernel, ncoeff)
    [respRaw, rGridRaw, sGridRaw] = rdRaw(data);
    rGridRawCal = rGridRaw - calRange;
    maskRange = rGridRawCal >= minRange & rGridRawCal <= maxRange;
    magnitudeRaw = abs(respRaw(maskRange, :));
    rangeGridFiltered = rGridRawCal(maskRange);
    
    mtiData = filter(mtiKernel, 1, data, [], 2);
    mtiData = mtiData(:, ncoeff:end);
    [respMTI, ~, sGridMTI] = rdFiltered(mtiData);
    
    magnitudeMTI_Full = abs(respMTI);
    [detMapFull, ~] = cfar_filter(magnitudeMTI_Full, 2, 5, 4, 20, 1e-2);
    
    detMapMasked = detMapFull(maskRange, :);
    magnitudeMTIMasked = magnitudeMTI_Full(maskRange, :);
end