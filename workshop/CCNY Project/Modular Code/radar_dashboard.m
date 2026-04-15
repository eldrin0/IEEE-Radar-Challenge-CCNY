function fns = radar_dashboard()
    fns.setup = @setup_ui;
    fns.update = @update_ui;
end

function handles = setup_ui()
    handles.fig = figure('Name', 'CCNY Radar Proximity Dashboard', 'Position', [50 50 1550 850], 'Color', 'w'); 
    
    handles.tlo = tiledlayout(2, 2, 'Padding', 'loose', 'TileSpacing', 'loose');

    handles.Raw = nexttile;
    handles.MTI = nexttile;
    handles.CFAR = nexttile;
    handles.Target = nexttile;

    % --- 1. Consistent Radar Look ---
    radarAxes = [handles.Raw, handles.MTI, handles.CFAR];
    for ax = radarAxes
        set(ax, 'Color', 'k', 'XColor', 'k', 'YColor', 'k');
        xlim(ax, [-15 15]); % Matching Raw Doppler axis
        ylim(ax, [0.2 10]); % Matching 10m range constraint
    end
    
    % --- 2. CFAR Grid (Dashed Lines) ---
    hold(handles.CFAR, 'on');
    grid(handles.CFAR, 'on');
    set(handles.CFAR, 'GridColor', [1 1 1], 'GridAlpha', 0.4, 'GridLineStyle', '--', 'Layer', 'top');

    % --- 3. Alert System Setup ---
    hold(handles.Target, 'on');
    handles.hSquare = rectangle('Parent', handles.Target, 'Position', [0.1, 0.25, 0.8, 0.55], 'Curvature', [0.1, 0.1], 'FaceColor', [0.8 0.8 0.8], 'LineWidth', 2);
    handles.hText = text(handles.Target, 0.5, 0.1, 'Searching...', 'HorizontalAlignment', 'center','FontSize', 14, 'FontWeight', 'bold');
    set(handles.Target, 'XTick', [], 'YTick', [], 'XColor', 'none', 'YColor', 'none');
    axis(handles.Target, [0 1 0 1]);
    
    colormap(handles.Raw, jet); colormap(handles.MTI, jet); colormap(handles.CFAR, parula);
end

function update_ui(handles, sGrid, rGrid, magRaw, magMTI, detMap, color, txt)
    %disp(handles.Target)
   

    % 1. Raw Response
    imagesc(handles.Raw, sGrid, rGrid, 10*log10(magRaw));
    axis(handles.Raw, 'xy');
    colorbar(handles.Raw);
    set(handles.Raw, 'CLim', [20 100]);
    title(handles.Raw, '1. Raw Response', 'FontSize', 14);

    % 2. MTI Filtered
    imagesc(handles.MTI, sGrid, rGrid, 10*log10(magMTI));
    axis(handles.MTI, 'xy');
    colorbar(handles.MTI);
    set(handles.MTI, 'CLim', [10 70]);
    title(handles.MTI, '2. MTI Filtered', 'FontSize', 14);

    % 3. CFAR Detections
    imagesc(handles.CFAR, sGrid, rGrid, double(detMap));
    axis(handles.CFAR, 'xy'); 
    colorbar(handles.CFAR); 
    set(handles.CFAR, 'CLim', [0 1]);
    title(handles.CFAR, '3. CFAR Detections', 'FontSize', 14);

    % 4. Alert Box
    set(handles.hSquare, 'FaceColor', color);
    set(handles.hText, 'String', txt);
    title(handles.Target, 'Target Proximity Status', 'FontSize', 14);
    
    drawnow limitrate;
end