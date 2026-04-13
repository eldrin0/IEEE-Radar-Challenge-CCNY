function fns = radar_dashboard()
    fns.setup = @setup_ui;
    fns.update = @update_ui;
end

function handles = setup_ui()
    handles.fig = figure('Name', 'CCNY Radar Proximity Dashboard', 'Position', [50 50 1550 850], 'Color', 'w'); 
    
    handles.tlo = tiledlayout(2, 2, 'Padding', 'loose', 'TileSpacing', 'loose');

    handles.axRaw = nexttile;
    handles.axMTI = nexttile;
    handles.axCFAR = nexttile;
    handles.axTarget = nexttile;

    % --- 1. Consistent Radar Look ---
    radarAxes = [handles.axRaw, handles.axMTI, handles.axCFAR];
    for ax = radarAxes
        set(ax, 'Color', 'k', 'XColor', 'k', 'YColor', 'k');
        xlim(ax, [-15 15]); % Matching Raw Doppler axis
        ylim(ax, [0.2 10]); % Matching 10m range constraint
    end

    % --- 2. Specific CFAR Grid (Dashed Lines) ---
    hold(handles.axCFAR, 'on');
    grid(handles.axCFAR, 'on');
    set(handles.axCFAR, 'GridColor', [1 1 1], 'GridAlpha', 0.4, 'GridLineStyle', '--', 'Layer', 'top');

    % --- 3. Alert System Setup ---
    hold(handles.axTarget, 'on');
    handles.hSquare = rectangle('Parent', handles.axTarget, 'Position', [0.1, 0.25, 0.8, 0.55], 'Curvature', [0.1, 0.1], 'FaceColor', [0.8 0.8 0.8], 'LineWidth', 2);
    handles.hText = text(handles.axTarget, 0.5, 0.1, 'Searching...', 'HorizontalAlignment', 'center','FontSize', 14, 'FontWeight', 'bold');
    set(handles.axTarget, 'XTick', [], 'YTick', [], 'XColor', 'none', 'YColor', 'none');
    axis(handles.axTarget, [0 1 0 1]);
    
    colormap(handles.axRaw, jet); colormap(handles.axMTI, jet); colormap(handles.axCFAR, parula);
end

function update_ui(handles, sGrid, rGrid, magRaw, magMTI, detMap, color, txt)
    % 1. Raw Response
    imagesc(handles.axRaw, sGrid, rGrid, 10*log10(magRaw));
    axis(handles.axRaw, 'xy');
    colorbar(handles.axRaw);
    set(handles.axRaw, 'CLim', [20 100]);
    title(handles.axRaw, '1. Raw Response', 'FontSize', 14);

    % 2. MTI Filtered
    imagesc(handles.axMTI, sGrid, rGrid, 10*log10(magMTI));
    axis(handles.axMTI, 'xy');
    colorbar(handles.axMTI);
    set(handles.axMTI, 'CLim', [10 70]);
    title(handles.axMTI, '2. MTI Filtered', 'FontSize', 14);

    % 3. CFAR Detections
    imagesc(handles.axCFAR, sGrid, rGrid, double(detMap));
    axis(handles.axCFAR, 'xy'); 
    colorbar(handles.axCFAR); 
    set(handles.axCFAR, 'CLim', [0 1]);
    title(handles.axCFAR, '3. CFAR Detections', 'FontSize', 14);

    % 4. Alert Box
    set(handles.hSquare, 'FaceColor', color);
    set(handles.hText, 'String', txt);
    title(handles.axTarget, 'Target Proximity Status', 'FontSize', 14);
    
    drawnow limitrate;
end