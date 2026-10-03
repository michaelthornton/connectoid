function fig = stim_plot_response_evidence(M, directions, label, opts)
% Where the two reported numbers come from, one panel per direction.
%
% Each dot is one trial's rate in the response window, in the organoid that was
% not stimulated. Dashed: the 95th percentile of that organoid's baseline
% windows. Solid with a band: the mean over trials and its 95% interval. Grey:
% the null's upper 97.5% bound, which the band must clear to count as
% connected.
arguments
    M (1,2) struct
    directions (1,2) string
    label (1,1) string
    opts.WindowSec (1,2) double = [0.030 0.100]
end
COL = [0 0 0; 0.62 0.11 0.11];
HIT = [0.10 0.45 0.85];

fig = figure('Units', 'pixels', 'Position', [60 60 1120 440], 'Color', 'w');
for p = 1:2
    m = M(p);
    r = m.trialRate(:)';
    n = numel(r);
    ax = axes(fig, 'Units', 'pixels', 'Position', [100 + (p - 1) * 540 96 420 240]);
    hold(ax, 'on');

    hi = r > m.thrSpS;
    plot(ax, find(~hi), r(~hi), 'o', 'MarkerSize', 4, 'Color', [0.62 0.62 0.62], ...
        'MarkerFaceColor', [0.82 0.82 0.82], 'LineWidth', 0.5);
    plot(ax, find(hi), r(hi), 'o', 'MarkerSize', 4.5, 'Color', HIT, ...
        'MarkerFaceColor', HIT, 'LineWidth', 0.5);

    patch(ax, [0 n + 1 n + 1 0], ...
        [m.respSpS - m.respCI, m.respSpS - m.respCI, ...
         m.respSpS + m.respCI, m.respSpS + m.respCI], COL(p, :), ...
        'FaceAlpha', 0.16, 'EdgeColor', 'none');
    plot(ax, [0 n + 1], [m.respSpS m.respSpS], '-', 'Color', COL(p, :), 'LineWidth', 1.6);
    plot(ax, [0 n + 1], [m.nullMeanP975 m.nullMeanP975], '-', ...
        'Color', [0.55 0.55 0.55], 'LineWidth', 1.2);
    plot(ax, [0 n + 1], [m.thrSpS m.thrSpS], '--', 'Color', HIT, 'LineWidth', 1.2);

    set(ax, 'XLim', [0 n + 1], 'FontSize', 11);
    stim_style(ax, 'FontSize', 11, 'LineWidth', 1.1);
    xlabel(ax, 'Trial', 'FontSize', 12);
    ylabel(ax, sprintf('Rate in the %g-%g ms window (spikes/s)', ...
        1000 * opts.WindowSec(1), 1000 * opts.WindowSec(2)), 'FontSize', 11);
    stim_left_title(ax, directions(p), 14);

    txt = sprintf(['%.0f%% of trials above the dashed line (p = %.3f)\n' ...
                   'window mean %+.0f%% of baseline, 95%% CI %+.0f to %+.0f%%\n' ...
                   '%s'], ...
        m.fracAbove, m.pFrac, m.deltaPct, m.deltaPct - m.deltaPctCI, ...
        m.deltaPct + m.deltaPctCI, callWord(m.ciSeparated));
    stim_norm_text(ax, 0.02, 1.34, txt, 'FontSize', 10, 'Color', [0.3 0.3 0.3], ...
        'VerticalAlignment', 'top');
end
annotation(fig, 'textbox', [0.07 0.005 0.9 0.045], 'EdgeColor', 'none', ...
    'FontSize', 10, 'Color', [0.35 0.35 0.35], 'VerticalAlignment', 'top', ...
    'String', sprintf(['Connectoid %s. Dashed: 95th percentile of the organoid''s own baseline windows. ' ...
                       'Solid with band: mean over trials, 95%% interval. Grey: the null''s upper 97.5%% bound.'], label));
end


function s = callWord(tf)
if tf
    s = 'connected: the interval clears the null';
else
    s = 'not connected: the interval overlaps the null';
end
end


function h = stim_left_title(ax, str, fontSize)
% A bold panel title flush with the left edge of the axes.
%
% Set Units and Position after creating: same ordering trap as stim_norm_text.
arguments
    ax
    str (1,1) string
    fontSize (1,1) double = 13.5
end
h = title(ax, str, 'FontSize', fontSize, 'FontWeight', 'bold', ...
          'HorizontalAlignment', 'left');
set(h, 'Units', 'normalized', 'Position', [0 1.02 0]);
end
