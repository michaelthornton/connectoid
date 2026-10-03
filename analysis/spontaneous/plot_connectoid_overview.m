function figures = plot_connectoid_overview(results, maxCorrelograms)
% The six QC figures for a two-organoid scan, in launcher save order.
%
%   1 electrode firing filter      4 accepted count correlograms
%   2 DBSCAN units and organoids   5 organoid activity
%   3 directed connections         6 inter-organoid synchrony
%
% results is the analyze_monosynaptic_connections struct with .activity and
% .synchrony attached. Nothing here computes a statistic.
arguments
    results struct
    maxCorrelograms (1,1) double {mustBePositive,mustBeInteger} = 12
end
figures = gobjects(0);
figures(end+1) = plotElectrodeFilter(results);
figures(end+1) = plotUnitsAndOrganoids(results);
figures(end+1) = plotDirectedConnections(results);
figures(end+1) = plotConnectionCounts(results, maxCorrelograms);
figures(end+1) = plotOrganoidActivity(results.activity);
figures(end+1) = plotOrganoidSynchrony(results.unitData, results.activity, results.synchrony);
end

function fig = plotElectrodeFilter(results)
E = results.electrodeData;
fig = figure('Color', 'w', 'Name', 'Step 1 - electrode firing filter', ...
    'WindowStyle', 'docked');
ax = axes(fig);
hold(ax, 'on');
scatter(ax, E.x(~E.passesFiringFilter), E.y(~E.passesFiringFilter), ...
    10, [0.82 0.82 0.82], 'filled', 'DisplayName', 'Excluded');
scatter(ax, E.x(E.passesFiringFilter), E.y(E.passesFiringFilter), ...
    18, [0.10 0.35 0.70], 'filled', 'DisplayName', 'Included');
formatArrayAxes(ax);
legend(ax, 'Location', 'best');
title(ax, sprintf('Electrode firing filter: %.3g to %.3g Hz', ...
    results.parameters.MinElectrodeFiringRateHz, ...
    results.parameters.MaxElectrodeFiringRateHz));
end

function fig = plotUnitsAndOrganoids(results)
E = results.electrodeData;
U = results.unitData;
fig = figure('Color', 'w', 'Name', 'Step 2 - DBSCAN units and organoids', ...
    'WindowStyle', 'docked');
ax = axes(fig);
hold(ax, 'on');
scatter(ax, E.x, E.y, 7, [0.88 0.88 0.88], 'filled');

organoidColors = lines(max(height(results.organoidRois), 1));
for iUnit = 1:height(U)
    memberIdx = U.memberMapIndices{iUnit};
    color = colorForOrganoid(U.organoidLabel(iUnit), results.organoidRois, organoidColors);
    scatter(ax, E.x(memberIdx), E.y(memberIdx), 34, color, 'filled', ...
        'MarkerEdgeColor', 'k');
    plot(ax, U.centroidX(iUnit), U.centroidY(iUnit), 'ko', ...
        'MarkerFaceColor', 'w', 'MarkerSize', 5);
    text(ax, U.centroidX(iUnit) + 10, U.centroidY(iUnit), ...
        string(U.unitId(iUnit)), 'FontSize', 8);
end
drawRois(ax, results.organoidRois, organoidColors);
formatArrayAxes(ax);
title(ax, sprintf('DBSCAN units: radius %.1f micrometers, min electrodes %d', ...
    results.parameters.DbscanRadiusUm, results.parameters.DbscanMinElectrodes));
subtitle(ax, sprintf('%d units have > %d spikes and are eligible for pair testing', ...
    nnz(U.passesPairSpikeFilter & U.organoidLabel ~= "unassigned"), ...
    results.parameters.MinUnitSpikesForPair));
end

function fig = plotDirectedConnections(results)
U = results.unitData;
C = results.connectionTable;
fig = figure('Color', 'w', 'Name', 'Step 3 - directed connections', ...
    'WindowStyle', 'docked');
ax = axes(fig);
hold(ax, 'on');
organoidColors = lines(max(height(results.organoidRois), 1));

for iUnit = 1:height(U)
    color = colorForOrganoid(U.organoidLabel(iUnit), results.organoidRois, organoidColors);
    scatter(ax, U.centroidX(iUnit), U.centroidY(iUnit), 65, color, ...
        'filled', 'MarkerEdgeColor', 'k');
    text(ax, U.centroidX(iUnit) + 10, U.centroidY(iUnit), ...
        string(U.unitId(iUnit)), 'FontSize', 8);
end

if ~isempty(C)
    % Arrow weight is relative peak prominence. The one-count floor only
    % keeps zero-baseline pairs finite; acceptance is unaffected.
    arrowProminence = C.peakCount ./ max(C.baselineMeanCount, 1);
    maxArrowProminence = max(arrowProminence);
    for iConnection = 1:height(C)
        source = C.sourceUnit(iConnection);
        target = C.targetUnit(iConnection);
        x1 = U.centroidX(source);
        y1 = U.centroidY(source);
        dx = U.centroidX(target) - x1;
        dy = U.centroidY(target) - y1;
        width = 0.75 + 2.25 * ...
            arrowProminence(iConnection) / maxArrowProminence;
        arrowColor = connectionColor(C.sourceOrganoid(iConnection), ...
            C.targetOrganoid(iConnection), results.organoidRois, organoidColors);
        quiver(ax, x1, y1, dx, dy, 0, 'Color', arrowColor, ...
            'LineWidth', width, 'MaxHeadSize', 0.25);
    end
else
    text(ax, 0.5, 0.5, 'No pairs passed the count and peak/baseline rules', ...
        'Units', 'normalized', 'HorizontalAlignment', 'center');
end

drawRois(ax, results.organoidRois, organoidColors);
formatArrayAxes(ax);
title(ax, sprintf(['Directed +%.1f to +%.1f ms peaks; ', ...
    'peak >= %d counts and >= %.1f x mean off-peak baseline'], ...
    results.parameters.PeakLatencyRangeMs(1), ...
    results.parameters.PeakLatencyRangeMs(2), ...
    results.parameters.MinPeakCount, ...
    results.parameters.PeakToBaselineMultiplier), 'Interpreter', 'none');
subtitle(ax, ...
    'Arrow width: peak count / max(mean off-peak baseline, 1 count)', ...
    'Interpreter', 'none');
end

function fig = plotConnectionCounts(results, maxPlots)
C = results.connectionTable;
fig = figure('Color', 'w', 'Name', 'Step 4 - connection counts', ...
    'WindowStyle', 'docked');

if isempty(C)
    ax = axes(fig);
    axis(ax, 'off');
    text(ax, 0.5, 0.5, 'No pairs passed the count and peak/baseline rules', ...
        'HorizontalAlignment', 'center');
    return;
end

[~, order] = sort(C.peakCount, 'descend');
order = order(1:min(maxPlots, numel(order)));
nPlots = numel(order);
nColumns = min(3, nPlots);
nRows = ceil(nPlots / nColumns);
layout = tiledlayout(fig, nRows, nColumns, 'TileSpacing', 'compact', ...
    'Padding', 'compact');

for iPlot = 1:nPlots
    row = C(order(iPlot), :);
    mask = results.pairCountTable.sourceUnit == row.sourceUnit & ...
        results.pairCountTable.targetUnit == row.targetUnit;
    counts = results.pairCountTable(mask, :);
    ax = nexttile(layout);
    bar(ax, counts.latencyMs, counts.count, 1, ...
        'FaceColor', [0.15 0.35 0.65], 'EdgeColor', 'none');
    hold(ax, 'on');
    yline(ax, row.baselineMeanCount * results.parameters.PeakToBaselineMultiplier, ...
        'r--');
    yline(ax, results.parameters.MinPeakCount, 'm--');
    xline(ax, 0, '-', 'Color', [0.65 0.65 0.65]);
    xline(ax, results.parameters.PeakLatencyRangeMs(1), 'g--');
    xline(ax, results.parameters.PeakLatencyRangeMs(2), 'g--');
    xline(ax, row.peakLatencyMs, 'k:');
    box(ax, 'off');
    xlabel(ax, 'Target relative to source (ms)');
    ylabel(ax, 'Pair count');
    title(ax, sprintf('Unit %d to %d | peak %d at %.1f ms, baseline %.1f', ...
        row.sourceUnit, row.targetUnit, row.peakCount, ...
        row.peakLatencyMs, row.baselineMeanCount), 'Interpreter', 'none');
end
title(layout, 'Accepted peaked short-latency count correlograms');
subtitle(layout, sprintf(['Green: +%.1f to +%.1f ms peak search; ', ...
    'red: %.1f x mean off-peak baseline; magenta: minimum peak count'], ...
    results.parameters.PeakLatencyRangeMs(1), ...
    results.parameters.PeakLatencyRangeMs(2), ...
    results.parameters.PeakToBaselineMultiplier));
end

function fig = plotOrganoidActivity(activity)
% Compare organoid metrics and their mean smoothed firing-rate traces.
fig = figure('Color', 'w', 'Name', 'Step 5 - organoid activity', ...
    'WindowStyle', 'docked');
layout = tiledlayout(fig, 2, 4, 'TileSpacing', 'compact', ...
    'Padding', 'compact');
colors = lines(max(numel(activity.organoidLabels), 2));
S = activity.summaryTable;

metricVariables = {'meanUnitFiringRateHz', 'peakFrequencyHz', ...
    'meanPeakFiringRateHz', 'peakToMeanRateRatio'};
metricTitles = {'Mean unit firing rate', 'Detected-peak frequency', ...
    'Mean detected-peak firing rate', 'Peak-to-mean rate ratio'};
metricYLabels = {'Hz/unit', 'peaks/s', 'Hz/unit', 'peak / mean'};
for iMetric = 1:numel(metricVariables)
    ax = nexttile(layout);
    values = S.(metricVariables{iMetric});
    bars = bar(ax, values, 'FaceColor', 'flat');
    bars.CData = colors(1:height(S), :);
    ax.XTick = 1:height(S);
    ax.XTickLabel = cellstr(S.organoidLabel);
    ylabel(ax, metricYLabels{iMetric});
    title(ax, metricTitles{iMetric});
    ymax = max(values, [], 'omitnan');
    ylim(ax, [0 max(1, 1.2 * ymax)]);
    box(ax, 'off');
end

ax = nexttile(layout, [1 4]);
hold(ax, 'on');
yMax = max(activity.organoidMeanRateHz, [], 'all');
if isempty(yMax) || ~isfinite(yMax) || yMax <= 0
    yMax = 1;
end

meanHandles = gobjects(numel(activity.organoidLabels), 1);
meanLabels = strings(numel(activity.organoidLabels), 1);
for iOrganoid = 1:numel(activity.organoidLabels)
    meanHandles(iOrganoid) = plot(ax, activity.traceTimeSec, ...
        activity.organoidMeanRateHz(:, iOrganoid), ...
        'Color', colors(iOrganoid, :), 'LineWidth', 1.8);
    meanLabels(iOrganoid) = activity.organoidLabels(iOrganoid) + " mean";
end

peakHandle = gobjects(0);
for iOrganoid = 1:numel(activity.organoidLabels)
    h = plot(ax, activity.peakTimesSec{iOrganoid}, activity.peakRatesHz{iOrganoid}, ...
        'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 5, 'LineStyle', 'none');
    if isempty(peakHandle) && ~isempty(activity.peakTimesSec{iOrganoid})
        peakHandle = h;
    else
        h.HandleVisibility = 'off';
    end
end

xlim(ax, [0 activity.parameters.RecordingDurationSec]);
ylim(ax, [0 1.05 * yMax]);
xlabel(ax, 'Time (s)');
ylabel(ax, 'Smoothed firing rate (Hz/unit)');
thresholdText = strings(numel(activity.organoidLabels), 1);
for iOrganoid = 1:numel(activity.organoidLabels)
    thresholdText(iOrganoid) = sprintf('%s: %d peaks, mean %.2f Hz/unit', ...
        activity.organoidLabels(iOrganoid), S.peakCount(iOrganoid), ...
        S.meanPeakFiringRateHz(iOrganoid));
end
title(ax, 'Organoid mean firing rates and detected peaks');
subtitle(ax, sprintf('Minimum prominence %.2f Hz; minimum distance %.1f s | %s', ...
    activity.parameters.PeakMinimumProminenceHz, ...
    activity.parameters.PeakMinimumDistanceSec, strjoin(thresholdText, ' | ')), ...
    'Interpreter', 'none');
box(ax, 'off');
if isempty(peakHandle)
    legend(ax, meanHandles, cellstr(meanLabels), 'Location', 'eastoutside');
else
    legend(ax, [meanHandles; peakHandle], ...
        cellstr([meanLabels; "Detected peaks"]), 'Location', 'eastoutside');
end
title(layout, 'Unit-level activity by organoid');
end

function fig = plotOrganoidSynchrony(unitData, activity, synchrony)
% Rasters, detected events, and directional delay histograms.
fig = figure('Color', 'w', 'Name', 'Step 6 - inter-organoid synchrony', ...
    'WindowStyle', 'docked');
layout = tiledlayout(fig, 3, 2, 'TileSpacing', 'compact', ...
    'Padding', 'compact');
labels = synchrony.organoidLabels;
colors = lines(2);

axRaster = nexttile(layout, [1 2]);
hold(axRaster, 'on');
rasterHandles = gobjects(2, 1);
groupCenters = nan(2, 1);
row = 0;
for iOrganoid = 1:2
    rows = find(string(unitData.organoidLabel) == labels(iOrganoid));
    firstRow = row + 1;
    for iUnit = 1:numel(rows)
        row = row + 1;
        spikeTimes = double(unitData.spikeTimesSec{rows(iUnit)}(:));
        plot(axRaster, spikeTimes, repmat(row, numel(spikeTimes), 1), '.', ...
            'Color', colors(iOrganoid, :), 'MarkerSize', 3, ...
            'HandleVisibility', 'off');
    end
    if row >= firstRow
        groupCenters(iOrganoid) = mean(firstRow:row);
    end
    rasterHandles(iOrganoid) = plot(axRaster, nan, nan, '.', ...
        'Color', colors(iOrganoid, :), 'MarkerSize', 10);
    if iOrganoid == 1 && row > 0
        yline(axRaster, row + 0.5, 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8, ...
            'HandleVisibility', 'off');
    end
end
xlim(axRaster, [0 activity.parameters.RecordingDurationSec]);
ylim(axRaster, [0 max(row + 1, 1)]);
yticks(axRaster, groupCenters(isfinite(groupCenters)));
yticklabels(axRaster, cellstr(labels(isfinite(groupCenters))));
ylabel(axRaster, 'Clustered units');
title(axRaster, 'Combined spike raster');
legend(axRaster, rasterHandles, cellstr(labels), 'Location', 'eastoutside');
box(axRaster, 'off');

axRate = nexttile(layout, [1 2]);
hold(axRate, 'on');
rateHandles = gobjects(2, 1);
for iOrganoid = 1:2
    rateHandles(iOrganoid) = plot(axRate, activity.traceTimeSec, ...
        activity.organoidMeanRateHz(:, iOrganoid), 'Color', colors(iOrganoid, :), ...
        'LineWidth', 1.4);
    if isfield(activity, 'timingPeakTimesSec')
        eventTimesSec = activity.timingPeakTimesSec{iOrganoid};
        eventRatesHz = interp1(activity.traceTimeSec, ...
            activity.organoidMeanRateHz(:, iOrganoid), eventTimesSec, 'linear');
    else
        eventTimesSec = activity.peakTimesSec{iOrganoid};
        eventRatesHz = activity.peakRatesHz{iOrganoid};
    end
    plot(axRate, eventTimesSec, eventRatesHz, ...
        'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 4, 'LineStyle', 'none', ...
        'HandleVisibility', 'off');
end
xlim(axRate, [0 activity.parameters.RecordingDurationSec]);
xlabel(axRate, 'Time (s)');
ylabel(axRate, 'Smoothed firing rate (Hz/unit)');
title(axRate, 'Mean firing-rate traces and timing-refined synchrony events');
legend(axRate, rateHandles, cellstr(labels), 'Location', 'eastoutside');
box(axRate, 'off');

delayHistogram = synchrony.postTriggerHistogramTable;
D = synchrony.directionalSummaryTable;
for iDirection = 1:height(D)
    axDelay = nexttile(layout);
    rows = delayHistogram.sourceOrganoid == D.sourceOrganoid(iDirection) & ...
        delayHistogram.targetOrganoid == D.targetOrganoid(iDirection);
    histogram = delayHistogram(rows, :);
    bar(axDelay, histogram.delayBinCenterSec, histogram.peakPairCount, 1, ...
        'FaceColor', colors(iDirection, :), 'EdgeColor', 'none');
    hold(axDelay, 'on');
    yline(axDelay, D.nullGlobalMaximum95(iDirection), '--', ...
        '95% circular-shift maximum', 'Color', [0.85 0.2 0.2], ...
        'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off');
    xlim(axDelay, [0 synchrony.parameters.MaxLagSec]);
    xlabel(axDelay, sprintf('Delay to %s peak after source (s)', D.targetOrganoid(iDirection)));
    ylabel(axDelay, 'All target-peak pairs');
    title(axDelay, sprintf(['%s to %s: %d pairs at %.3g s; ', ...
        'global-shift p = %.3g'], ...
        D.sourceOrganoid(iDirection), D.targetOrganoid(iDirection), ...
        D.modalPairCount(iDirection), D.modalDelaySec(iDirection), ...
        D.empiricalPValue(iDirection)), ...
        'Interpreter', 'none');
    box(axDelay, 'off');
end
title(layout, 'Inter-organoid all-pairs event timing with circular-shift null');
end

function color = colorForOrganoid(label, organoidRois, organoidColors)
idx = find(organoidRois.label == label, 1);
if isempty(idx)
    color = [0.45 0.45 0.45];
else
    color = organoidColors(idx, :);
end
end

function color = connectionColor(sourceLabel, targetLabel, organoidRois, organoidColors)
if sourceLabel == targetLabel
    color = [0.35 0.35 0.35];
else
    color = colorForOrganoid(sourceLabel, organoidRois, organoidColors);
end
end

function drawRois(ax, organoidRois, organoidColors)
for iOrganoid = 1:height(organoidRois)
    position = normalizedRectanglePosition(organoidRois.position(iOrganoid, :));
    rectangle(ax, 'Position', position, 'EdgeColor', organoidColors(iOrganoid, :), ...
        'LineWidth', 2, 'LineStyle', '--');
    text(ax, position(1), position(2) + position(4), organoidRois.label(iOrganoid), ...
        'Color', organoidColors(iOrganoid, :), 'FontWeight', 'bold', ...
        'Interpreter', 'none');
end

function position = normalizedRectanglePosition(position)
bounds = [min(position(1), position(1) + position(3)), ...
    max(position(1), position(1) + position(3)), ...
    min(position(2), position(2) + position(4)), ...
    max(position(2), position(2) + position(4))];
position = [bounds(1), bounds(3), bounds(2) - bounds(1), bounds(4) - bounds(3)];
end
end

function formatArrayAxes(ax)
axis(ax, 'equal');
axis(ax, 'ij');
axis(ax, 'tight');
box(ax, 'off');
xlabel(ax, 'x (\mum)');
ylabel(ax, 'y (\mum)');
end
