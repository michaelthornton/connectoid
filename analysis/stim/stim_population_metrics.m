function M = stim_population_metrics(S, organoid, opts)
% Response of one organoid in the window after each train, against a null
% drawn from its own stimulus-free baseline.
%
%   fracAbove  % of trials whose window rate exceeds the 95th percentile of the
%              baseline windows. Its no-effect value is 5 by construction.
%   deltaPct   mean over trials as a % change from the baseline mean.
%              ciSeparated is true when its 95% interval clears the central 95%
%              of the null's.
%
% The null is every window position on a 1 ms grid across the baseline,
% resampled at the trial count.
arguments
    S struct
    organoid struct
    opts.WindowSec (1,2) double = [0.030 0.100]
    opts.Iterations (1,1) double {mustBeInteger,mustBePositive} = 20000
    opts.GridStepSec (1,1) double {mustBePositive} = 0.001
    opts.Seed (1,1) double = 1
end
w = opts.WindowSec;
assert(w(2) > w(1), 'stim:Window', 'WindowSec must be increasing.');
span = w(2) - w(1);
t0 = S.onsets(1);                       % stimulus-free baseline is [0, t0)
nTrials = numel(S.onsets);

ts = sort(S.t(ismember(S.el, organoid.electrodes)));

nGrid = ceil((t0 - w(2)) / opts.GridStepSec);
assert(nGrid > 1000, 'stim:ShortBaseline', ...
    'Only %d baseline windows fit before the first train.', nGrid);
pool = sort(stim_grid_window_counts(ts, nGrid, opts.GridStepSec, w) / span);

% Order statistic rather than an interpolated percentile, so the threshold is
% implementation independent.
thr = pool(ceil(0.95 * numel(pool)));
base = mean(pool);

resp = stim_trial_window_counts(ts, S.onsets, w) / span;
frac = 100 * mean(resp > thr);
delta = 100 * (mean(resp) - base) / base;

rng(opts.Seed, 'twister');
draws = pool(randi(numel(pool), opts.Iterations, nTrials));
nullFrac = 100 * mean(draws > thr, 2);
nullDelta = 100 * (mean(draws, 2) - base) / base;
pval = @(obs, nul) (1 + sum(nul >= obs)) / (opts.Iterations + 1);

nullMeans = mean(draws, 2);
loObs = mean(resp) - 1.96 * std(resp) / sqrt(nTrials);
hiNull = quantile(nullMeans, 0.975);
% False-positive rate of each rule, measured on the null itself.
fprDelta = mean(nullDelta > quantile(nullDelta, 0.95));
loNull = mean(draws, 2) - 1.96 * std(draws, 0, 2) / sqrt(nTrials);
fprCI = mean(loNull > hiNull);
nullP95Pct  = 100 * (quantile(nullMeans, 0.95) - base) / base;
nullP975Pct = 100 * (hiNull - base) / base;

M = struct('respondingOrganoid', organoid.type, ...
           'nElectrodes', organoid.nElectrodes, 'nTrials', nTrials, ...
           'nNullWindows', numel(pool), ...
           'baselineSpS', base, 'thrSpS', thr, 'respSpS', mean(resp), ...
           'fracAbove', frac, 'deltaPct', delta, 'trialRate', resp, ...
           'pFrac', pval(frac, nullFrac), 'pDelta', pval(delta, nullDelta), ...
           'ciSeparated', loObs > hiNull, 'fprDelta', fprDelta, ...
           'fprCI', fprCI, ...
           'nullP95Pct', nullP95Pct, 'nullP975Pct', nullP975Pct, ...
           'nullMeanP975', hiNull, ...
           'respCI', 1.96 * std(resp) / sqrt(nTrials), ...
           'deltaPctCI', 100 * 1.96 * std(resp) / sqrt(nTrials) / base);
end


function n = stim_trial_window_counts(ts, onsets, w)
% Spikes in [onset+w(1), onset+w(2)) for each onset.
%
% The trains are 1 s apart and the window is 70 ms, so the windows never
% overlap and one histcounts over interleaved edges does it: the odd bins are
% the windows. A trailing edge keeps the last window half-open; histcounts
% closes its final bin.
ts = sort(ts(:));
onsets = onsets(:);
e = reshape([onsets + w(1), onsets + w(2)]', [], 1);
assert(all(diff(e) > 0), 'stim:OverlappingWindows', ...
    'The response windows overlap, so the interleaved-edge trick is invalid here.');
c = histcounts(ts, [e; e(end) + 1]);
n = c(1:2:end - 1)';
end


function n = stim_grid_window_counts(ts, nGrid, step, w)
% Counts in the same window slid along a uniform grid from zero.
%
% Window i starts at (i-1)*step, i = 1..nGrid. The windows overlap (1 ms step,
% 70 ms window), so this differences a cumulative count rather than searching
% per edge. Exact because w is a whole number of grid steps.
k0 = round(w(1) / step);
k1 = round(w(2) / step);
assert(abs(k0 * step - w(1)) < 1e-12 && abs(k1 * step - w(2)) < 1e-12, ...
    'stim:GridStep', 'The window bounds must be whole multiples of the grid step.');
edges = (0:(nGrid + k1)) * step;
cnt = histcounts(ts, edges);
C = [0 cumsum(cnt)];                  % C(j) = spikes in [0, (j-1)*step)
i = (1:nGrid)';
n = C(i + k1)' - C(i + k0)';
end
