function stim_style(ax, opts)
% The house axis style: no box, thick spines, outward ticks.
arguments
    ax
    opts.FontSize (1,1) double = 13
    opts.LineWidth (1,1) double = 1.4
end
box(ax, 'off');
set(ax, 'TickDir', 'out', 'LineWidth', opts.LineWidth, 'FontSize', opts.FontSize, ...
        'Layer', 'top', 'TickLength', [0.012 0.012], 'Color', 'none');
end
