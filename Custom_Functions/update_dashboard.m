function [hImgRaw, hImgMTI, hImgCFAR] = update_dashboard(hImgRaw, hImgMTI, hImgCFAR, axRaw, axMTI, axCFAR, sGridRaw, rangeGridFiltered, magnitudeRaw, sGridMTI, magnitudeMTIMasked, detMapMasked, hSquare, storedColor, hText, storedText)
    if isempty(hImgRaw)
        % Plot 1: Raw
        hImgRaw = imagesc(axRaw, sGridRaw, rangeGridFiltered, 10*log10(magnitudeRaw));
        title(axRaw, '1. Raw Response'); colorbar(axRaw); axis(axRaw, 'xy');
        xlabel(axRaw, 'Speed (m/s)'); ylabel(axRaw, 'Range (m)');
        
        % Plot 2: MTI
        hImgMTI = imagesc(axMTI, sGridMTI, rangeGridFiltered, 10*log10(magnitudeMTIMasked));
        title(axMTI, '2. MTI Filtered'); colorbar(axMTI); axis(axMTI, 'xy');
        xlabel(axMTI, 'Speed (m/s)');
        
        % Plot 3: CFAR
        hImgCFAR = imagesc(axCFAR, sGridMTI, rangeGridFiltered, detMapMasked);
        title(axCFAR, '3. CFAR Detections'); colorbar(axCFAR); axis(axCFAR, 'xy');
        xlabel(axCFAR, 'Speed (m/s)');
    else
        % Update Data Only
        set(hImgRaw, 'CData', 10*log10(magnitudeRaw));
        set(hImgMTI, 'CData', 10*log10(magnitudeMTIMasked));
        set(hImgCFAR, 'CData', detMapMasked);
    end
    
    % Update Alert UI
    set(hSquare, 'FaceColor', storedColor);
    set(hText, 'String', storedText);
    
    drawnow limitrate;
end