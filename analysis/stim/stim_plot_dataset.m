function fig = stim_plot_dataset(recs, R, M, label, opts)
% One connectoid, both stimulation directions, on one sheet.
%
% A column per stimulation file and a row per organoid: trials as a heatmap,
% and beneath it the same trials as a population rate in % of that organoid's
% own baseline, on the same time axis.
%
% Heatmaps are scaled per panel and rates are per organoid, because the two
% organoids differ in routed electrode count. The two undriven rate panels
% share limits so the directions can be compared; the driven panels share a
% separate pair.
arguments
    recs (1,2) struct
    R table
    M (1,2) struct
    label (1,1) string
    opts.SpanSec (1,2) double = [-0.100 0.400]
    opts.BinSec (1,1) double = 0.020
    opts.CurveBinSec (1,1) double = 0.010
    opts.WindowSec (1,2) double = [0.030 0.100]
    opts.PreBlankSec (1,1) double = 0.005
    opts.ScalePct (1,1) double = 99
    opts.Gamma (1,1) double = 0.6
    opts.SmoothBins (1,1) double = 0.45
end
COL = [0 0 0; 0.62 0.11 0.11];      % row 1, row 2
MARK = [0.10 0.45 0.85];            % the blanked stimulation window
types = R.type(:)';
cm = stim_maxwell_colormap();

hedges = stim_anchored_edges(opts.SpanSec, opts.BinSec, opts.WindowSec(1));
cedges = stim_anchored_edges(opts.SpanSec, opts.CurveBinSec, opts.WindowSec(1));
hblank = stim_blank_mask(hedges, opts.PreBlankSec, opts.WindowSec(1));
cblank = stim_blank_mask(cedges, opts.PreBlankSec, opts.WindowSec(1));
tCurve = 500 * (cedges(1:end-1) + cedges(2:end));       % bin centres, ms
xlim0 = 1000 * opts.SpanSec;

% Everything is computed first so the shared limits see the whole sheet.
panel = struct('organoid', cell(2, 2));
for p = 1:2
    S = recs(p);
    [undriven, driven] = stim_organoids(S, R);
    for j = 1:2
        o = driven; if o.type ~= types(j), o = undriven; end
        ts = sort(S.t(ismember(S.el, o.electrodes)));
        counts = stim_trial_matrix(ts, S.onsets, hedges);
        rate = stim_trial_matrix(ts, S.onsets, cedges) / opts.CurveBinSec;
        t0 = S.onsets(1);
        baseHz = numel(ts(ts < t0)) / t0;
        assert(baseHz > 0, 'stim:NoBaseline', ...
            '%s: the %s organoid has no spikes before the first train.', ...
            S.recording, o.type);
        pct = 100 * rate / baseHz;
        mu = mean(pct, 1); sem = std(pct, 0, 1) / sqrt(size(pct, 1));
        mu(cblank) = NaN; sem(cblank) = NaN;
        panel(p, j).organoid = o;
        panel(p, j).isDriven = o.type == driven.type;
        panel(p, j).counts = counts;
        panel(p, j).mu = mu;
        panel(p, j).sem = sem;
        panel(p, j).stim = driven.type;
    end
end
dLim = bandLimit(panel, true);
uLim = bandLimit(panel, false);

FW = 1480; FH = 900;
X0 = 118; PW = 548; PGAP = 150; CBW = 15;
HH = 186; CH = 128; VGAP = 4; GRPGAP = 46;
yTop = FH - 96;
fig = figure('Units', 'pixels', 'Position', [40 40 FW FH], 'Color', 'w');

for p = 1:2
    S = recs(p);
    x = X0 + (p - 1) * (PW + PGAP);
    for j = 1:2
        q = panel(p, j);
        yh = yTop - (j - 1) * (HH + VGAP + CH + GRPGAP) - HH;
        yc = yh - VGAP - CH;

        % --- trials
        vmax = prctile(reshape(q.counts(:, ~hblank), [], 1), opts.ScalePct);
        if ~(vmax > 0), vmax = 1; end
        [Mf, blankFine] = stim_display_smooth(q.counts, hblank, opts.SmoothBins);
        up = size(Mf, 2) / numel(hblank);
        tFine = 1000 * (hedges(1) + opts.BinSec * ((0:size(Mf, 2) - 1) + 0.5) / up);
        % Gamma on the data, since MATLAB has no PowerNorm. max() drops NaN,
        % so the blank is re-masked after.
        Mg = vmax * (max(Mf, 0) / vmax) .^ opts.Gamma;
        Mg(:, blankFine) = NaN;
        ax = axes(fig, 'Units', 'pixels', 'Position', [x yh PW HH]);
        im = imagesc(ax, tFine, 1:size(Mg, 1), Mg);
        set(im, 'AlphaData', ~isnan(Mg));
        colormap(ax, cm); clim(ax, [0 vmax]);
        set(ax, 'YDir', 'reverse', 'XLim', xlim0, 'XTickLabel', [], ...
            'YTick', [1 size(Mg, 1)], 'FontSize', 10, 'Color', 'w', ...
            'Box', 'on', 'LineWidth', 1, 'TickDir', 'out');
        stim_norm_text(ax, -0.075, 0.5, q.organoid.type, 'Rotation', 90, ...
            'FontSize', 15, 'Color', COL(j, :), 'FontWeight', weightOf(q.isDriven), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
        stim_norm_text(ax, -0.118, 0.5, drivenWord(q.isDriven), 'Rotation', 90, ...
            'FontSize', 9.5, 'Color', [0.45 0.45 0.45], ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
        stim_norm_text(ax, 0.985, 0.06, sprintf('%d electrodes', q.organoid.nElectrodes), ...
            'FontSize', 9, 'Color', 'w', 'HorizontalAlignment', 'right', ...
            'VerticalAlignment', 'bottom');
        % Ticks go through the same gamma, so the labels stay in spikes/bin.
        tv = [0 vmax/4 vmax/2 vmax];
        cb = colorbar(ax, 'Units', 'pixels', 'Position', [x + PW + 12 yh CBW HH], ...
            'Ticks', vmax * (tv / vmax) .^ opts.Gamma, ...
            'TickLabels', compose('%.0f', tv), 'FontSize', 9, 'LineWidth', 0.8);
        if j == 1
            cb.Label.String = sprintf('spikes / %d ms bin', round(1000 * opts.BinSec));
            cb.Label.FontSize = 9;
        end

        % --- the same trials as a rate
        lim = uLim; if q.isDriven, lim = dLim; end
        ax = axes(fig, 'Units', 'pixels', 'Position', [x yc PW CH]);
        hold(ax, 'on');
        patch(ax, 1000 * opts.WindowSec([1 2 2 1]), lim([1 1 2 2]), MARK, ...
            'FaceAlpha', 0.10, 'EdgeColor', 'none');
        patch(ax, 1000 * [-opts.PreBlankSec opts.WindowSec(1) opts.WindowSec(1) -opts.PreBlankSec], ...
            lim([1 1 2 2]), [0.5 0.5 0.5], 'FaceAlpha', 0.12, 'EdgeColor', 'none');
        yline(ax, 100, ':', 'Color', [0.45 0.45 0.45], 'LineWidth', 0.9);
        drawBand(ax, tCurve, q.mu, q.sem, COL(j, :));
        set(ax, 'XLim', xlim0, 'YLim', lim, 'FontSize', 10);
        stim_style(ax, 'FontSize', 10, 'LineWidth', 1);
        ylabel(ax, {'Rate', '(% of baseline)'}, 'FontSize', 10);
        if j == 2
            xlabel(ax, 'Time from first pulse (ms)', 'FontSize', 12);
        else
            set(ax, 'XTickLabel', []);
        end
        if ~q.isDriven
            stim_norm_text(ax, 0.985, 0.90, statLine(M(p)), 'FontSize', 9.5, ...
                'Color', COL(j, :), 'HorizontalAlignment', 'right', ...
                'VerticalAlignment', 'top');
        end
    end
    header(fig, FW, x, PW, "Stimulate " + panel(p, 1).stim, ...
        sprintf('%s   %d trials', S.recording, numel(S.onsets)), ...
        COL(find(types == panel(p, 1).stim, 1), :));
end

annotation(fig, 'textbox', [X0 / FW 0.004 0.92 0.072], 'EdgeColor', 'none', ...
    'FontSize', 10, 'Color', [0.35 0.35 0.35], 'VerticalAlignment', 'top', ...
    'String', {sprintf(['Connectoid %s. Trials in acquisition order. Grey band: the %g to +%g ms blank, ' ...
                        'masked everywhere. Blue band: the %g to %g ms response window.'], ...
                       label, -1000 * opts.PreBlankSec, 1000 * opts.WindowSec(1), ...
                       1000 * opts.WindowSec(1), 1000 * opts.WindowSec(2)), ...
               sprintf(['Heatmaps are spikes per %d ms bin over each organoid''s routed electrodes, each scaled 0 to its own %dth percentile. ' ...
                        'Rate panels are the mean over trials, shaded +/- SEM, in %% of that organoid''s own pre-stimulus baseline.'], ...
                       round(1000 * opts.BinSec), opts.ScalePct)});
end


function drawBand(ax, x, mu, sem, col)
% Drawn in runs, so the blank stays a gap rather than a line across it.
good = ~isnan(mu);
d = diff([false, good, false]);
starts = find(d == 1);
stops = find(d == -1) - 1;
for s = 1:numel(starts)
    i = starts(s):stops(s);
    if numel(i) < 2, continue; end
    patch(ax, [x(i) fliplr(x(i))], [mu(i) - sem(i) fliplr(mu(i) + sem(i))], col, ...
        'FaceAlpha', 0.22, 'EdgeColor', 'none');
    plot(ax, x(i), mu(i), 'Color', col, 'LineWidth', 1.4);
end
end


function lim = bandLimit(panel, wantDriven)
% One pair of limits for the two panels of the same kind, rounded outward and
% not anchored at zero.
hi = -inf; lo = inf;
for p = 1:2
    for j = 1:2
        if panel(p, j).isDriven ~= wantDriven, continue; end
        hi = max(hi, max(panel(p, j).mu + panel(p, j).sem, [], 'omitnan'));
        lo = min(lo, min(panel(p, j).mu - panel(p, j).sem, [], 'omitnan'));
    end
end
step = 10 ^ floor(log10(max(hi - lo, 1)) - 0.3);
lim = [floor(lo / step) * step, ceil(hi / step) * step];
if ~(lim(2) > lim(1)), lim = [0 100]; end
end


function s = statLine(M)
% The two reported numbers for this direction.
mark = ""; if M.ciSeparated, mark = "  connected"; end
s = sprintf('%+.0f%% in window%s\n%.0f%% of trials respond', M.deltaPct, mark, M.fracAbove);
end


function w = weightOf(isDriven)
if isDriven, w = 'bold'; else, w = 'normal'; end
end


function s = drivenWord(isDriven)
if isDriven, s = 'driven'; else, s = 'undriven'; end
end


function header(fig, FW, x, PW, title1, title2, col)
annotation(fig, 'textbox', [(x + PW/2 - 160) / FW 0.955 320 / FW 0.040], ...
    'String', title1, 'EdgeColor', 'none', 'FontSize', 17, 'Color', col, ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
annotation(fig, 'textbox', [(x + PW/2 - 160) / FW 0.928 320 / FW 0.028], ...
    'String', title2, 'EdgeColor', 'none', 'FontSize', 10, ...
    'Color', [0.5 0.5 0.5], 'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle');
end


function edges = stim_anchored_edges(span, binSec, w0)
% Bin edges on the grid w0 + k*binSec, clipped to span, so w0 (the end of the
% blank) lands on an edge.
%
% A bin straddling the end of the blank would mix artifact into the first
% response bin. The far edge of the window is not on the grid, so that line
% falls inside a bin.
arguments
    span (1,2) double
    binSec (1,1) double {mustBePositive}
    w0 (1,1) double
end
assert(span(2) > span(1), 'stim:Span', 'span must be increasing.');
k0 = ceil((span(1) - w0) / binSec);
k1 = floor((span(2) - w0) / binSec);
edges = w0 + binSec * (k0:k1);
end


function m = stim_blank_mask(edges, preBlankSec, windowStartSec)
% Bins overlapping the blank, so they can be masked.
%
% Overlap, not bin centre: accumulated rounding puts a centre at
% -5.000000000000051, which fails a >= -5 test and leaves one artifact bin in.
m = (1000 * edges(2:end) > -1000 * preBlankSec) & ...
    (1000 * edges(1:end-1) < 1000 * windowStartSec);
end


function [Mf, blankFine] = stim_display_smooth(M, blank, sigmaBins, up)
% Upsample and smooth a trials-by-bins image along time only.
%
% Along time only; smoothing across trials would invent structure. Time is
% upsampled `up`-fold with a Gaussian of sigmaBins original bins, and
% sigmaBins = 0 returns the raw bins. The blank is held out of the kernel and
% re-applied as NaN, so no artifact crosses its edge.
arguments
    M double
    blank logical
    sigmaBins (1,1) double = 0.45
    up (1,1) double {mustBeInteger,mustBePositive} = 8
end
nb = size(M, 2);
assert(numel(blank) == nb, 'stim:BlankSize', ...
    'blank has %d entries for %d bins.', numel(blank), nb);
if sigmaBins <= 0 || up <= 1
    Mf = M; blankFine = blank(:)'; return
end
x = (1:nb) - 0.5;                                  % original bin centres
xf = (0:nb * up - 1) / up + 0.5 / up;              % fine grid centres
blankFine = repelem(blank(:)', up);
W = exp(-0.5 * ((xf - x(:)) / sigmaBins) .^ 2);     % nb-by-(nb*up)
W(blank(:), :) = 0;                                 % the blank contributes nothing
nrm = sum(W, 1);
nrm(nrm == 0) = 1;
Mf = (M * W) ./ nrm;
Mf(:, blankFine) = NaN;
end


function cmap = stim_maxwell_colormap(n)
% The MaxWell activity-map colour scale: black at zero, through blue, white in
% the middle, to red at the top.
%
% Sampled from the MaxWell firing-rate maps, so panels drawn here match the
% ones their software exports.
arguments
    n (1,1) double {mustBePositive,mustBeInteger} = 256
end
stops = [  2   4   6;  26  57 101;  53 106 173;  85 154 215; 110 182 235
         132 194 238; 157 207 241; 183 220 245; 205 231 248; 232 244 252
         254 251 247; 250 228 201; 247 205 155; 244 182 111; 241 159  71
         238 128  49; 234  56  36] / 255;
cmap = interp1(linspace(0, 1, size(stops, 1)), stops, linspace(0, 1, n)');
end


function M = stim_trial_matrix(ts, onsets, edges)
% trials-by-bins spike counts, edges relative to each onset.
%
% A loop, not interleaved edges: the span runs -200 to +950 ms while the trains
% are 1 s apart, so consecutive trials' spans overlap and the edges would not
% be monotonic.
ts = sort(ts(:));
M = zeros(numel(onsets), numel(edges) - 1);
for k = 1:numel(onsets)
    lo = onsets(k) + edges(1);
    hi = onsets(k) + edges(end);
    seg = ts(ts >= lo & ts < hi) - onsets(k);
    M(k, :) = histcounts(seg, edges);
end
end
