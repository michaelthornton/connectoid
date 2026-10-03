function activity = analyze_organoid_activity(unitData, varargin)
% Unit-rate summaries and sparse MFR peaks.
% Clustered unit spike times only; no voltage data is reopened.

p = inputParser;
p.addRequired('unitData', @istable);
p.addParameter('RecordingDurationSec', [], @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('OrganoidLabels', {'Organoid 1', 'Organoid 2'}, ...
    @(x) iscellstr(x) || isstring(x));
p.addParameter('RateBinSec', 0.01, @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('RateGaussianSigmaSec', 0.5, ...
    @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('TimingGaussianSigmaSec', 0.05, ...
    @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('TimingRefinementWindowSec', 0.5, ...
    @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('PeakMinimumDistanceSec', 2, ...
    @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('PeakMinimumProminenceHz', 1, ...
    @(x) isnumeric(x) && isscalar(x) && x > 0);
p.parse(unitData, varargin{:});
opts = p.Results;

requiredVariables = {'unitId', 'spikeTimesSec', 'organoidLabel'};
missingVariables = setdiff(requiredVariables, unitData.Properties.VariableNames);
if ~isempty(missingVariables)
    error('analyze_organoid_activity:InvalidUnitTable', ...
        'unitData is missing: %s.', strjoin(missingVariables, ', '));
end
if isempty(which('mxw.networkActivity.computeNetworkAct'))
    error('analyze_organoid_activity:MissingToolbox', ...
        'Add the MaxWell MATLAB toolbox before running this function.');
end

labels = string(opts.OrganoidLabels(:));
if isempty(labels) || any(strlength(labels) == 0) || numel(unique(labels)) ~= numel(labels)
    error('analyze_organoid_activity:InvalidOrganoidLabels', ...
        'Organoid labels must be nonempty and unique.');
end

assigned = ismember(string(unitData.organoidLabel), labels);
unitData = unitData(assigned, :);
traceTimeSec = (0:opts.RateBinSec:opts.RecordingDurationSec).';
unitRateHz = zeros(numel(traceTimeSec), height(unitData));
unitMeanFiringRateHz = zeros(height(unitData), 1);

for iUnit = 1:height(unitData)
    spikeTimes = sort(double(unitData.spikeTimesSec{iUnit}(:)));
    spikeTimes = spikeTimes(spikeTimes >= 0 & spikeTimes <= opts.RecordingDurationSec);
    unitMeanFiringRateHz(iUnit) = numel(spikeTimes) / opts.RecordingDurationSec;
    spikeStruct = struct('time', spikeTimes, ...
        'channel', repmat(double(unitData.unitId(iUnit)), numel(spikeTimes), 1));
    rate = mxw.networkActivity.computeNetworkAct(spikeStruct, ...
        'BinSize', opts.RateBinSec, ...
        'GaussianSigma', opts.RateGaussianSigmaSec, ...
        'MinValue', 0, 'MaxValue', opts.RecordingDurationSec);
    unitRateHz(:, iUnit) = rate.firingRate(:);
end

nOrganoids = numel(labels);
organoidMeanRateHz = zeros(numel(traceTimeSec), nOrganoids);
timingRateHz = zeros(numel(traceTimeSec), nOrganoids);
nUnits = zeros(nOrganoids, 1);
meanUnitFiringRateHz = zeros(nOrganoids, 1);
peakCount = zeros(nOrganoids, 1);
peakFrequencyHz = zeros(nOrganoids, 1);
meanPeakFiringRateHz = zeros(nOrganoids, 1);
peakToMeanRateRatio = zeros(nOrganoids, 1);
peakTimesSec = cell(nOrganoids, 1);
peakRatesHz = cell(nOrganoids, 1);
peakProminencesHz = cell(nOrganoids, 1);
timingPeakTimesSec = cell(nOrganoids, 1);
timingPeakRatesHz = cell(nOrganoids, 1);
peakTables = cell(nOrganoids, 1);

for iOrganoid = 1:nOrganoids
    inOrganoid = string(unitData.organoidLabel) == labels(iOrganoid);
    nUnits(iOrganoid) = nnz(inOrganoid);
    if nUnits(iOrganoid) == 0
        peakTimesSec{iOrganoid} = zeros(0, 1);
        peakRatesHz{iOrganoid} = zeros(0, 1);
        peakProminencesHz{iOrganoid} = zeros(0, 1);
        timingPeakTimesSec{iOrganoid} = zeros(0, 1);
        timingPeakRatesHz{iOrganoid} = zeros(0, 1);
        peakTables{iOrganoid} = emptyPeakTable();
        continue;
    end

    organoidMeanRateHz(:, iOrganoid) = mean(unitRateHz(:, inOrganoid), 2);
    meanUnitFiringRateHz(iOrganoid) = mean(unitMeanFiringRateHz(inOrganoid));
    [peakRatesHz{iOrganoid}, peakTimesSec{iOrganoid}, ~, peakProminencesHz{iOrganoid}] = findpeaks( ...
        organoidMeanRateHz(:, iOrganoid), traceTimeSec, ...
        'MinPeakDistance', opts.PeakMinimumDistanceSec, ...
        'MinPeakProminence', opts.PeakMinimumProminenceHz);
    allSpikeTimesSec = vertcat(unitData.spikeTimesSec{inOrganoid});
    allSpikeTimesSec = sort(double(allSpikeTimesSec(:)));
    allSpikeTimesSec = allSpikeTimesSec(allSpikeTimesSec >= 0 & ...
        allSpikeTimesSec <= opts.RecordingDurationSec);
    populationSpikes = struct('time', allSpikeTimesSec, ...
        'channel', ones(numel(allSpikeTimesSec), 1));
    timingRate = mxw.networkActivity.computeNetworkAct(populationSpikes, ...
        'BinSize', opts.RateBinSec, ...
        'GaussianSigma', opts.TimingGaussianSigmaSec, ...
        'MinValue', 0, 'MaxValue', opts.RecordingDurationSec);
    timingRateHz(:, iOrganoid) = timingRate.firingRate(:) / nUnits(iOrganoid);
    [timingPeakTimesSec{iOrganoid}, timingPeakRatesHz{iOrganoid}] = ...
        refinePeakTimes(peakTimesSec{iOrganoid}, traceTimeSec, ...
        timingRateHz(:, iOrganoid), opts.TimingRefinementWindowSec);
    peakCount(iOrganoid) = numel(peakRatesHz{iOrganoid});
    peakFrequencyHz(iOrganoid) = peakCount(iOrganoid) / opts.RecordingDurationSec;
    if peakCount(iOrganoid) > 0
        meanPeakFiringRateHz(iOrganoid) = mean(peakRatesHz{iOrganoid});
    end
    if meanUnitFiringRateHz(iOrganoid) > 0
        peakToMeanRateRatio(iOrganoid) = meanPeakFiringRateHz(iOrganoid) / ...
            meanUnitFiringRateHz(iOrganoid);
    else
        peakToMeanRateRatio(iOrganoid) = NaN;
    end
    peakTables{iOrganoid} = table( ...
        repmat(labels(iOrganoid), peakCount(iOrganoid), 1), ...
        (1:peakCount(iOrganoid)).', peakTimesSec{iOrganoid}, ...
        timingPeakTimesSec{iOrganoid}, peakRatesHz{iOrganoid}, ...
        timingPeakRatesHz{iOrganoid}, ...
        peakProminencesHz{iOrganoid}, repmat(opts.PeakMinimumProminenceHz, ...
        peakCount(iOrganoid), 1), ...
        'VariableNames', {'organoidLabel', 'peakNumber', 'peakTimeSec', ...
        'timingPeakTimeSec', 'peakFiringRateHz', 'timingPeakFiringRateHz', ...
        'peakProminenceHz', 'minimumPeakProminenceHz'});
end

summaryTable = table(labels, nUnits, meanUnitFiringRateHz, peakCount, ...
    peakFrequencyHz, meanPeakFiringRateHz, peakToMeanRateRatio, ...
    'VariableNames', {'organoidLabel', 'nUnits', 'meanUnitFiringRateHz', ...
    'peakCount', 'peakFrequencyHz', 'meanPeakFiringRateHz', ...
    'peakToMeanRateRatio'});
unitSummaryTable = table(unitData.unitId, string(unitData.organoidLabel), ...
    unitMeanFiringRateHz, 'VariableNames', {'unitId', 'organoidLabel', ...
    'meanFiringRateHz'});

activity = struct();
activity.parameters = rmfield(opts, 'unitData');
activity.organoidLabels = labels;
activity.traceTimeSec = traceTimeSec;
activity.unitIds = unitData.unitId;
activity.unitOrganoidLabels = string(unitData.organoidLabel);
activity.unitRateHz = unitRateHz;
activity.organoidMeanRateHz = organoidMeanRateHz;
activity.timingRateHz = timingRateHz;
activity.summaryTable = summaryTable;
activity.unitSummaryTable = unitSummaryTable;
activity.peakTimesSec = peakTimesSec;
activity.peakRatesHz = peakRatesHz;
activity.peakProminencesHz = peakProminencesHz;
activity.timingPeakTimesSec = timingPeakTimesSec;
activity.timingPeakRatesHz = timingPeakRatesHz;
activity.peakTable = vertcat(peakTables{:});
end

function [refinedTimesSec, refinedRatesHz] = refinePeakTimes( ...
        candidateTimesSec, traceTimeSec, rateHz, refinementWindowSec)
refinedTimesSec = zeros(numel(candidateTimesSec), 1);
refinedRatesHz = zeros(numel(candidateTimesSec), 1);
for iPeak = 1:numel(candidateTimesSec)
    inWindow = abs(traceTimeSec - candidateTimesSec(iPeak)) <= refinementWindowSec;
    [refinedRatesHz(iPeak), localIndex] = max(rateHz(inWindow));
    windowTimesSec = traceTimeSec(inWindow);
    refinedTimesSec(iPeak) = windowTimesSec(localIndex);
end
end

function peakTable = emptyPeakTable()
peakTable = table(strings(0, 1), zeros(0, 1), zeros(0, 1), zeros(0, 1), ...
    zeros(0, 1), zeros(0, 1), zeros(0, 1), zeros(0, 1), ...
    'VariableNames', {'organoidLabel', 'peakNumber', 'peakTimeSec', ...
    'timingPeakTimeSec', 'peakFiringRateHz', 'timingPeakFiringRateHz', ...
    'peakProminenceHz', 'minimumPeakProminenceHz'});
end
