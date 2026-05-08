function [hMainFig, axRaw, axMTI, axCFAR, axTarget, hSquare, hText] = setup_ui()
    hMainFig = figure('Name', 'Radar Dashboard', 'Position', [100 100 1200 800], 'Color', 'w');
    tlo = tiledlayout(2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
    
    % Prepare axes
    axRaw = nexttile; axMTI = nexttile; axCFAR = nexttile; axTarget = nexttile;
    
    % Tile 4: Diagnostic Alert System
    hold(axTarget, 'on');
    hSquare = rectangle('Parent', axTarget, 'Position', [0.1, 0.2, 0.8, 0.6], ...
        'Curvature', [0.1, 0.1], 'FaceColor', [0.8 0.8 0.8], 'LineWidth', 2);
    hText = text(axTarget, 0.5, 0.1, 'Searching...', 'HorizontalAlignment', 'center', ...
        'FontSize', 12, 'FontWeight', 'bold');
    title(axTarget, 'Target Proximity Status');
    set(axTarget, 'XTick', [], 'YTick', [], 'XColor', 'none', 'YColor', 'none');
    axis(axTarget, [0 1 0 1]);
end