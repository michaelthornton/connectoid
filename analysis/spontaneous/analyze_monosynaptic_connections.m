function results = analyze_monosynaptic_connections(recordingInput, varargin)
% Count short-latency directed spike pairs.
%
%   results = analyze_monosynaptic_connections(dataFile, Name, Value, ...)
%
% Loads a MaxWell H5 with mxw.fileManager, filters electrodes by firing rate,
% clusters them with DBSCAN, and cross-correlates with pxcorr from the
% toolbox's bundled UltraMegaSort. Positive latency is target after source.
%
% Counts are kept over +/-LatencyRangeMs (default -11 to +11 ms); the peak is
% searched only in PeakLatencyRangeMs (default +1 to +11 ms) and the baseline
% is the mean of the bins outside it (-11 to -1 ms). A pair passes at
% MinPeakCount and PeakToBaselineMultiplier times that baseline. No shuffle
% correction, no z-scoring. Units outside the ROIs, or at or below
% MinUnitSpikesForPair spikes, stay in the QC table but are not tested.
%
% Duplicate units are rejected by testing for duplication directly: shared
% member electrodes, or centroids inside DuplicateRadiusUm, together with a
% correlogram dominated by the |lag| < 1 ms band. MaxPeakEfficacy is a backstop
% rather than the duplicate test: peak efficacy scales with coupling strength
% and not with duplication, so a genuine local connection is strong. The
% launchers set it to 0.5, where a pair whose target follows more than half the
% source's spikes is rejected whatever the spatial test says. peakEfficacy and
% nearZeroLagCount are reported either way.
%
% For tests, recordingInput may be an electrode table with electrode, x, y
% and spikeTimesSec; RecordingDurationSec must then be set.

p = inputParser;
p.FunctionName = mfilename;
p.addRequired('recordingInput', @(x) ischar(x) || isstring(x) || istable(x));
p.addParameter('WellID', 1, @(x) isscalar(x) && isnumeric(x) && x >= 1);
p.addParameter('RecordingNumber', 1, @(x) isscalar(x) && isnumeric(x) && x >= 1);
p.addParameter('RecordingDurationSec', [], @(x) isempty(x) || (isscalar(x) && isnumeric(x) && x > 0));
p.addParameter('MinElectrodeFiringRateHz', 0.01, @(x) isscalar(x) && isnumeric(x) && x >= 0);
p.addParameter('MaxElectrodeFiringRateHz', Inf, @(x) isscalar(x) && isnumeric(x) && x > 0);
p.addParameter('DbscanRadiusUm', 60, @(x) isscalar(x) && isnumeric(x) && x > 0);
p.addParameter('DbscanMinElectrodes', 2, @(x) isscalar(x) && isnumeric(x) && x >= 1 && mod(x, 1) == 0);
p.addParameter('MergeWindowMs', 1.5, @(x) isscalar(x) && isnumeric(x) && x >= 0);
p.addParameter('MinUnitSpikesForPair', 50, @(x) isscalar(x) && isnumeric(x) && ...
    isfinite(x) && x >= 0 && mod(x, 1) == 0);
p.addParameter('CorrelationBinMs', 0.5, @(x) isscalar(x) && isnumeric(x) && isfinite(x) && x > 0);
p.addParameter('LatencyRangeMs', [1 11], @(x) isnumeric(x) && numel(x) == 2 && ...
    all(isfinite(x)) && x(1) > 0 && x(2) >= x(1));
p.addParameter('PeakLatencyRangeMs', [1 11], @(x) isnumeric(x) && numel(x) == 2 && ...
    all(isfinite(x)) && x(1) > 0 && x(2) >= x(1));
p.addParameter('MinPeakCount', 5, @(x) isscalar(x) && isnumeric(x) && ...
    isfinite(x) && x >= 1 && mod(x, 1) == 0);
p.addParameter('PeakToBaselineMultiplier', 5, @(x) isscalar(x) && isnumeric(x) && x > 0);
p.addParameter('MaxPeakEfficacy', Inf, @(x) isscalar(x) && isnumeric(x) && x > 0);
% Two units are the same neuron if they share electrodes or sit inside the
% clustering radius. Default matches DbscanRadiusUm.
p.addParameter('DuplicateRadiusUm', 60, @(x) isscalar(x) && isnumeric(x) && x >= 0);
% And a correlogram whose |lag| < 1 ms band, per bin, exceeds this multiple of
% the peak bin is co-firing rather than transmission.
p.addParameter('MaxNearZeroLagRatio', 2, @(x) isscalar(x) && isnumeric(x) && x > 0);
p.addParameter('OrganoidLabels', {'Organoid 1', 'Organoid 2'}, @(x) iscellstr(x) || isstring(x));
p.addParameter('AllowStimulationEvents', false, @(x) islogical(x) && isscalar(x));
% What to call each organoid on the ROI picker, when the labels are positional.
p.addParameter('OrganoidDescriptions', string.empty, @(x) iscellstr(x) || isstring(x));
p.addParameter('OrganoidRoiPositions', [], @(x) isempty(x) || ...
    (isnumeric(x) && size(x, 2) == 4 && all(isfinite(x), 'all')));
p.addParameter('InteractiveRoiSelection', true, @(x) islogical(x) && isscalar(x));
p.addParameter('UseSquareRois', true, @(x) islogical(x) && isscalar(x));
p.addParameter('AssignAllUnitsToSingleOrganoid', false, ...
    @(x) islogical(x) && isscalar(x));
p.parse(recordingInput, varargin{:});
opts = p.Results;

if exist('dbscan', 'file') ~= 2
    error('analyze_monosynaptic_connections:MissingDbscan', ...
        'MATLAB''s dbscan function is required.');
end
if exist('pxcorr', 'file') ~= 2
    error('analyze_monosynaptic_connections:MissingPxcorr', ...
        ['pxcorr was not found. Add the MaxWell matlab folder recursively ', ...
         'before calling this function.']);
end

stageTimer = tic;
fprintf('Reading spikes...');
[electrodeData, recordingDurationSec, sourceMetadata] = ...
    loadElectrodeData(recordingInput, opts);
fprintf(' %.1f s\n', toc(stageTimer));
opts.RecordingDurationSec = recordingDurationSec;

if opts.MaxElectrodeFiringRateHz < opts.MinElectrodeFiringRateHz
    error('analyze_monosynaptic_connections:InvalidFiringRateRange', ...
        'MaxElectrodeFiringRateHz must be at least MinElectrodeFiringRateHz.');
end

electrodeData.nSpikes = cellfun(@numel, electrodeData.spikeTimesSec);
electrodeData.firingRateHz = electrodeData.nSpikes ./ recordingDurationSec;
electrodeData.passesFiringFilter = ...
    electrodeData.firingRateHz >= opts.MinElectrodeFiringRateHz & ...
    electrodeData.firingRateHz <= opts.MaxElectrodeFiringRateHz;

stageTimer = tic;
fprintf('Clustering %d electrodes past the %.3g Hz filter...', ...
    sum(electrodeData.passesFiringFilter), opts.MinElectrodeFiringRateHz);
unitData = clusterUnits(electrodeData, opts);
fprintf(' %d units, %.1f s\n', height(unitData), toc(stageTimer));
organoidRois = resolveOrganoidRois(electrodeData, opts);
unitData.organoidLabel = assignUnitsToOrganoids(unitData, organoidRois);
unitData.passesPairSpikeFilter = unitData.nSpikes > opts.MinUnitSpikesForPair;

stageTimer = tic;
[pairCountTable, pairSummaryTable] = computeDirectedPairCounts(unitData, opts);
fprintf(' %.1f s\n', toc(stageTimer));
connectionTable = pairSummaryTable(pairSummaryTable.isConnection, :);
connectionSummaryTable = summarizeConnectionsByOrganoid( ...
    connectionTable, organoidRois.label);

parameters = rmfield(opts, {'recordingInput', 'RecordingDurationSec'});
parameters.RecordingDurationSec = recordingDurationSec;

provenance = struct();
provenance.AnalysisFunction = mfilename;
provenance.AnalysisVersion = '1.3.0';
provenance.CreatedAt = char(datetime('now', 'TimeZone', 'local', ...
    'Format', 'yyyy-MM-dd HH:mm:ss Z'));
provenance.MatlabVersion = version;
provenance.MatlabRelease = version('-release');
provenance.Source = sourceMetadata;
provenance.CorrelationFunction = which('pxcorr');
provenance.ClusteringFunction = which('dbscan');

results = struct();
results.parameters = parameters;
results.provenance = provenance;
results.organoidRois = organoidRois;
results.electrodeData = electrodeData;
results.unitData = unitData;
results.pairCountTable = pairCountTable;
results.pairSummaryTable = pairSummaryTable;
results.connectionTable = connectionTable;
results.connectionSummaryTable = connectionSummaryTable;
end

function summaryTable = summarizeConnectionsByOrganoid(connectionTable, labels)
% One label yields a within-organoid count; two labels yield the paired summary.
if isscalar(labels)
    sourceOrganoid = labels(1);
    targetOrganoid = labels(1);
    relationship = labels(1) + " within";
elseif numel(labels) == 2
    sourceOrganoid = [labels(1); labels(2); labels(1); labels(2)];
    targetOrganoid = [labels(1); labels(2); labels(2); labels(1)];
    relationship = [labels(1) + " within"; labels(2) + " within"; ...
        labels(1) + " to " + labels(2); labels(2) + " to " + labels(1)];
else
    error('analyze_monosynaptic_connections:InvalidOrganoidLabels', ...
        'Connectivity summaries support one or two organoid labels.');
end
nConnections = zeros(numel(relationship), 1);
for iRelationship = 1:numel(relationship)
    nConnections(iRelationship) = nnz( ...
        connectionTable.sourceOrganoid == sourceOrganoid(iRelationship) & ...
        connectionTable.targetOrganoid == targetOrganoid(iRelationship));
end
summaryTable = table(relationship, sourceOrganoid, targetOrganoid, nConnections);
end

function [electrodeData, durationSec, metadata] = loadElectrodeData(recordingInput, opts)
if istable(recordingInput)
    requiredVariables = {'electrode', 'x', 'y', 'spikeTimesSec'};
    missingVariables = setdiff(requiredVariables, recordingInput.Properties.VariableNames);
    if ~isempty(missingVariables)
        error('analyze_monosynaptic_connections:InvalidElectrodeTable', ...
            'Electrode table is missing: %s.', strjoin(missingVariables, ', '));
    end
    if isempty(opts.RecordingDurationSec)
        error('analyze_monosynaptic_connections:MissingDuration', ...
            'RecordingDurationSec is required when recordingInput is a table.');
    end
    electrodeData = recordingInput(:, requiredVariables);
    electrodeData.electrode = double(electrodeData.electrode(:));
    electrodeData.x = double(electrodeData.x(:));
    electrodeData.y = double(electrodeData.y(:));
    electrodeData.spikeTimesSec = cellfun(@(x) sort(double(x(:))), ...
        electrodeData.spikeTimesSec, 'UniformOutput', false);
    durationSec = double(opts.RecordingDurationSec);
    metadata = struct('Type', 'synthetic_or_preloaded_table', 'DataFile', '');
    return;
end

dataFile = char(recordingInput);
if ~isfile(dataFile)
    error('analyze_monosynaptic_connections:MissingDataFile', ...
        'Data file does not exist: %s', dataFile);
end

nStimEvents = countStimulationEvents(dataFile, opts.WellID);
if nStimEvents > 0 && ~opts.AllowStimulationEvents
    error('analyze_monosynaptic_connections:StimulationRecording', ...
        ['%s holds %d stimulation events, so it is a stimulation recording ' ...
         'rather than a network scan. Point this analysis at the network ' ...
         'scan, or set AllowStimulationEvents to true to run it anyway.'], ...
        dataFile, nStimEvents);
end

recording = mxw.fileManager(dataFile, opts.WellID);
if opts.RecordingNumber > recording.nRecordings
    error('analyze_monosynaptic_connections:InvalidRecordingNumber', ...
        'RecordingNumber %d exceeds the %d recordings in the file.', ...
        opts.RecordingNumber, recording.nRecordings);
end

rec = recording.fileObj(opts.RecordingNumber);
fs = double(rec.samplingFreq);
durationSec = double(rec.dataLenSamples) / fs;

spikeFrames = normalizeSpikeFramesByElectrode( ...
    recording.extractedSpikes(opts.RecordingNumber).frameno);
electrode = double(recording.processedMap.electrode(:));
x = double(recording.processedMap.xpos(:));
y = double(recording.processedMap.ypos(:));

if ~iscell(spikeFrames) || numel(spikeFrames) ~= numel(electrode)
    error('analyze_monosynaptic_connections:UnexpectedMaxwellData', ...
        'MaxWell extractedSpikes does not align with processedMap.');
end

firstFrame = resolveFirstFrame(rec.firstFrameNum, spikeFrames);
spikeTimesSec = cellfun(@(frames) cleanSpikeTimes(frames, firstFrame, fs, durationSec), ...
    spikeFrames(:), 'UniformOutput', false);
electrodeData = table(electrode, x, y, spikeTimesSec);

fileInfo = dir(dataFile);
metadata = struct();
metadata.Type = 'maxwell_h5';
metadata.DataFile = dataFile;
metadata.FileBytes = fileInfo.bytes;
metadata.FileModified = fileInfo.date;
metadata.WellID = opts.WellID;
metadata.RecordingNumber = opts.RecordingNumber;
metadata.SamplingFrequencyHz = fs;
metadata.FirstFrame = firstFrame;
end

function nEvents = countStimulationEvents(dataFile, wellID)
% Stimulation events in the file, read from the dataset size alone.
nEvents = 0;
try
    info = h5info(dataFile, sprintf('/data_store/data%04d/events', wellID - 1));
catch
    return;
end
if ~isempty(info.Dataspace.Size)
    nEvents = double(info.Dataspace.Size(1));
end
end

function spikeFrames = normalizeSpikeFramesByElectrode(spikeFrames)
if iscell(spikeFrames)
    spikeFrames = spikeFrames(:);
elseif isstruct(spikeFrames) && isfield(spikeFrames, 'frameno')
    spikeFrames = {spikeFrames.frameno}.';
elseif isnumeric(spikeFrames)
    spikeFrames = {spikeFrames(:)};
else
    error('analyze_monosynaptic_connections:UnexpectedSpikeFormat', ...
        'Unsupported MaxWell extractedSpikes frameno format.');
end
end

function firstFrame = resolveFirstFrame(firstFrameValue, spikeFrames)
if isnumeric(firstFrameValue) && isscalar(firstFrameValue) && isfinite(firstFrameValue)
    firstFrame = double(firstFrameValue);
    return;
end

firstFrame = Inf;
for iElectrode = 1:numel(spikeFrames)
    frames = double(spikeFrames{iElectrode}(:));
    frames = frames(isfinite(frames));
    if ~isempty(frames)
        firstFrame = min(firstFrame, min(frames));
    end
end
if isinf(firstFrame)
    firstFrame = 0;
end
warning('analyze_monosynaptic_connections:MissingFirstFrame', ...
    ['MaxWell did not provide a numeric first frame. Relative times are ', ...
     'measured from the first detected spike, matching the toolbox fallback.']);
end

function timesSec = cleanSpikeTimes(frames, firstFrame, fs, durationSec)
timesSec = (double(frames(:)) - firstFrame) ./ fs;
timesSec = sort(timesSec(isfinite(timesSec) & timesSec >= 0 & timesSec <= durationSec));
end

function unitData = clusterUnits(electrodeData, opts)
selectedIdx = find(electrodeData.passesFiringFilter);

if isempty(selectedIdx)
    unitData = emptyUnitTable();
    return;
end

clusterId = dbscan([electrodeData.x(selectedIdx), electrodeData.y(selectedIdx)], ...
    opts.DbscanRadiusUm, opts.DbscanMinElectrodes);
groupIds = unique(clusterId(clusterId > 0), 'stable');
nUnits = numel(groupIds);

unitId = (1:nUnits).';
memberMapIndices = cell(nUnits, 1);
memberElectrodes = cell(nUnits, 1);
spikeTimesSec = cell(nUnits, 1);
centroidX = zeros(nUnits, 1);
centroidY = zeros(nUnits, 1);
nSpikes = zeros(nUnits, 1);
firingRateHz = zeros(nUnits, 1);

for iUnit = 1:nUnits
    mapIdx = selectedIdx(clusterId == groupIds(iUnit));
    combinedTimes = sort(vertcat(electrodeData.spikeTimesSec{mapIdx}));
    mergedTimes = mergeNearbyDetections(combinedTimes, opts.MergeWindowMs / 1000);

    memberMapIndices{iUnit} = mapIdx(:).';
    memberElectrodes{iUnit} = electrodeData.electrode(mapIdx).';
    spikeTimesSec{iUnit} = mergedTimes;
    centroidX(iUnit) = median(electrodeData.x(mapIdx));
    centroidY(iUnit) = median(electrodeData.y(mapIdx));
    nSpikes(iUnit) = numel(mergedTimes);
    firingRateHz(iUnit) = nSpikes(iUnit) / opts.RecordingDurationSec;
end

unitData = table(unitId, memberMapIndices, memberElectrodes, spikeTimesSec, ...
    centroidX, centroidY, nSpikes, firingRateHz, ...
    'VariableNames', {'unitId', 'memberMapIndices', 'memberElectrodes', ...
    'spikeTimesSec', 'centroidX', 'centroidY', 'nSpikes', 'firingRateHz'});
end

function unitData = emptyUnitTable()
unitData = table(zeros(0,1), cell(0,1), cell(0,1), cell(0,1), ...
    zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), ...
    'VariableNames', {'unitId', 'memberMapIndices', 'memberElectrodes', ...
    'spikeTimesSec', 'centroidX', 'centroidY', 'nSpikes', 'firingRateHz'});
end

function mergedTimes = mergeNearbyDetections(timesSec, mergeWindowSec)
if isempty(timesSec)
    mergedTimes = zeros(0, 1);
    return;
end
timesSec = sort(timesSec(:));
mergedTimes = timesSec([true; diff(timesSec) > mergeWindowSec]);
end

function organoidRois = resolveOrganoidRois(electrodeData, opts)
labels = string(opts.OrganoidLabels(:));
positions = double(opts.OrganoidRoiPositions);

if isempty(labels) || any(strlength(labels) == 0) || ...
        numel(unique(labels)) ~= numel(labels) || any(labels == "unassigned")
    error('analyze_monosynaptic_connections:InvalidOrganoidLabels', ...
        'Organoid labels must be nonempty, unique, and cannot be "unassigned".');
end

described = string(opts.OrganoidDescriptions(:));
if isempty(described)
    described = labels;
elseif numel(described) ~= numel(labels)
    error('analyze_monosynaptic_connections:InvalidOrganoidDescriptions', ...
        'OrganoidDescriptions must have one entry per OrganoidLabels entry.');
end

if opts.AssignAllUnitsToSingleOrganoid
    if numel(labels) ~= 1
        error('analyze_monosynaptic_connections:InvalidSingleOrganoidAssignment', ...
            'AssignAllUnitsToSingleOrganoid requires exactly one OrganoidLabels entry.');
    end
    positions = arrayBounds(electrodeData.x, electrodeData.y);
end

if isempty(positions)
    if ~opts.InteractiveRoiSelection
        error('analyze_monosynaptic_connections:MissingOrganoidRois', ...
            'Provide OrganoidRoiPositions or enable InteractiveRoiSelection.');
    end
    positions = select_organoid_rois(electrodeData, labels, opts.UseSquareRois, ...
        described);
end

if size(positions, 1) ~= numel(labels)
    error('analyze_monosynaptic_connections:InvalidOrganoidRois', ...
        'OrganoidRoiPositions must have one row per OrganoidLabels entry.');
end
if any(positions(:, 3) == 0 | positions(:, 4) == 0)
    error('analyze_monosynaptic_connections:InvalidOrganoidRois', ...
        'Every organoid ROI must have nonzero width and height.');
end

organoidRois = table(labels, positions, described, ...
    'VariableNames', {'label', 'position', 'type'});
end

function position = arrayBounds(x, y)
% Display-only box over the whole array: one organoid needs no selection.
if isempty(x) || isempty(y)
    position = [0 0 1 1];
    return;
end
xMin = min(x);
xMax = max(x);
yMin = min(y);
yMax = max(y);
position = [xMin, yMin, max(xMax - xMin, eps), max(yMax - yMin, eps)];
end

function labels = assignUnitsToOrganoids(unitData, organoidRois)
labels = repmat("unassigned", height(unitData), 1);
membershipCount = zeros(height(unitData), 1);

for iOrganoid = 1:height(organoidRois)
    bounds = rectangleBounds(organoidRois.position(iOrganoid, :));
    inside = unitData.centroidX >= bounds(1) & unitData.centroidX <= bounds(2) & ...
        unitData.centroidY >= bounds(3) & unitData.centroidY <= bounds(4);
    membershipCount = membershipCount + inside;
    labels(inside & labels == "unassigned") = organoidRois.label(iOrganoid);
end

if any(membershipCount > 1)
    error('analyze_monosynaptic_connections:OverlappingOrganoidRois', ...
        '%d units fall inside overlapping organoid ROIs.', nnz(membershipCount > 1));
end
end

function bounds = rectangleBounds(position)
x1 = position(1);
x2 = position(1) + position(3);
y1 = position(2);
y2 = position(2) + position(4);
bounds = [min(x1, x2), max(x1, x2), min(y1, y2), max(y1, y2)];
end

function [pairCounts, pairSummary] = computeDirectedPairCounts(unitData, opts)
assignedRows = find(unitData.organoidLabel ~= "unassigned" & ...
    unitData.passesPairSpikeFilter);
nUnits = numel(assignedRows);
latencyRangeMs = double(opts.LatencyRangeMs(:).');
peakLatencyRangeMs = double(opts.PeakLatencyRangeMs(:).');
if peakLatencyRangeMs(1) < latencyRangeMs(1) || ...
        peakLatencyRangeMs(2) > latencyRangeMs(2)
    error('analyze_monosynaptic_connections:InvalidPeakLatencyRange', ...
        'PeakLatencyRangeMs must lie inside LatencyRangeMs.');
end
correlationFs = 1000 / opts.CorrelationBinMs;
maxLagSec = latencyRangeMs(2) / 1000;
maxLagBins = ceil(maxLagSec * correlationFs);
availableLatenciesMs = 1000 * ((-maxLagBins:maxLagBins) / correlationFs);
expectedLatenciesMs = availableLatenciesMs( ...
    availableLatenciesMs >= -latencyRangeMs(2) & ...
    availableLatenciesMs <= latencyRangeMs(2));
nLatencyBins = numel(expectedLatenciesMs);
if nLatencyBins == 0
    error('analyze_monosynaptic_connections:NoLatencyBins', ...
        'CorrelationBinMs does not place a bin inside LatencyRangeMs.');
end
peakSearchMask = expectedLatenciesMs >= peakLatencyRangeMs(1) & ...
    expectedLatenciesMs <= peakLatencyRangeMs(2);
baselineMask = expectedLatenciesMs <= -latencyRangeMs(1) | ...
    expectedLatenciesMs > peakLatencyRangeMs(2);
if ~any(peakSearchMask) || ~any(baselineMask)
    error('analyze_monosynaptic_connections:InsufficientLagBins', ...
        'The configured lag ranges do not provide both peak and baseline bins.');
end
peakSearchIndices = find(peakSearchMask);
% Lags inside +/-LatencyRangeMs(1) are neither peak nor baseline: counts
% there are the co-detection signature, so report them.
nearZeroMask = abs(expectedLatenciesMs) < latencyRangeMs(1);
nPairs = nUnits * max(nUnits - 1, 0);

sourceUnit = zeros(nPairs * nLatencyBins, 1);
targetUnit = zeros(nPairs * nLatencyBins, 1);
sourceOrganoid = strings(nPairs * nLatencyBins, 1);
targetOrganoid = strings(nPairs * nLatencyBins, 1);
latencyMs = zeros(nPairs * nLatencyBins, 1);
count = zeros(nPairs * nLatencyBins, 1);

summarySource = zeros(nPairs, 1);
summaryTarget = zeros(nPairs, 1);
summarySourceOrganoid = strings(nPairs, 1);
summaryTargetOrganoid = strings(nPairs, 1);
summarySourceNSpikes = zeros(nPairs, 1);
summaryTargetNSpikes = zeros(nPairs, 1);
peakLatencyMs = nan(nPairs, 1);
peakCount = zeros(nPairs, 1);
baselineMeanCount = zeros(nPairs, 1);
peakToBaselineRatio = nan(nPairs, 1);
peakEfficacy = nan(nPairs, 1);
nearZeroLagCount = zeros(nPairs, 1);
sharedElectrodes = zeros(nPairs, 1);
centroidDistanceUm = nan(nPairs, 1);
nearZeroLagRatio = nan(nPairs, 1);
passesEfficacyCap = false(nPairs, 1);
passesDuplicateTest = false(nPairs, 1);
passesNearZeroLag = false(nPairs, 1);
passesMinPeakCount = false(nPairs, 1);
passesPeakToBaseline = false(nPairs, 1);
isConnection = false(nPairs, 1);

pairIdx = 0;
countRow = 0;
fprintf('Correlating %d units, %d directed pairs', nUnits, nUnits * (nUnits - 1));
progressEvery = max(1, round(nUnits / 20));
for iSource = 1:nUnits
    if mod(iSource, progressEvery) == 0
        fprintf('.');
    end
    sourceRow = assignedRows(iSource);
    for iTarget = 1:nUnits
        targetRow = assignedRows(iTarget);
        if sourceRow == targetRow
            continue;
        end
        pairIdx = pairIdx + 1;
        sourceTimes = unitData.spikeTimesSec{sourceRow};
        targetTimes = unitData.spikeTimesSec{targetRow};

        thisCounts = zeros(1, nLatencyBins);
        if numel(sourceTimes) > 1 && numel(targetTimes) > 1
            % pxcorr follows xcorr ordering: target first, source second
            % makes positive lag mean target after source.
            [correlationHz, pxcorrLagsSec] = pxcorr( ...
                targetTimes, sourceTimes, correlationFs, maxLagSec);
            rawCounts = round(correlationHz * numel(targetTimes) / correlationFs);
            pxcorrLagsMs = 1000 * pxcorrLagsSec;
            for iLatency = 1:nLatencyBins
                [distance, nearestIdx] = min(abs(pxcorrLagsMs - expectedLatenciesMs(iLatency)));
                if distance <= max(eps(expectedLatenciesMs(iLatency)) * 8, 1e-9)
                    thisCounts(iLatency) = rawCounts(nearestIdx);
                end
            end
        end

        [thisPeak, relativePeakIdx] = max(thisCounts(peakSearchMask));
        peakIdx = peakSearchIndices(relativePeakIdx);
        thisBaseline = mean(thisCounts(baselineMask));
        if thisBaseline == 0
            if thisPeak > 0
                thisRatio = Inf;
            else
                thisRatio = NaN;
            end
        else
            thisRatio = thisPeak / thisBaseline;
        end
        % Fraction of source spikes followed by a peak-bin target spike,
        % above chance. A synapse is a few percent, a split unit far more.
        thisEfficacy = max(thisPeak - thisBaseline, 0) / max(numel(sourceTimes), 1);
        passesCount = thisPeak >= opts.MinPeakCount;
        passesRatio = thisPeak > 0 && ...
            thisPeak >= opts.PeakToBaselineMultiplier * thisBaseline;
        passesEfficacy = thisEfficacy <= opts.MaxPeakEfficacy;
        % Are these two units the same neuron? Shared electrodes or centroids
        % inside the clustering radius, plus co-firing at |lag| < 1 ms.
        thisShared = numel(intersect(unitData.memberElectrodes{sourceRow}, ...
                                     unitData.memberElectrodes{targetRow}));
        thisDistance = hypot(unitData.centroidX(sourceRow) - unitData.centroidX(targetRow), ...
                             unitData.centroidY(sourceRow) - unitData.centroidY(targetRow));
        thisNearZero = sum(thisCounts(nearZeroMask));
        nNearZeroBins = max(sum(nearZeroMask), 1);
        if thisPeak > 0
            thisNearZeroRatio = (thisNearZero / nNearZeroBins) / thisPeak;
        else
            thisNearZeroRatio = Inf;
        end
        passesDuplicate = thisShared == 0 && thisDistance > opts.DuplicateRadiusUm;
        passesNearZero = thisNearZeroRatio < opts.MaxNearZeroLagRatio;
        passes = passesCount && passesRatio && passesEfficacy && ...
            passesDuplicate && passesNearZero;

        summarySource(pairIdx) = unitData.unitId(sourceRow);
        summaryTarget(pairIdx) = unitData.unitId(targetRow);
        summarySourceOrganoid(pairIdx) = unitData.organoidLabel(sourceRow);
        summaryTargetOrganoid(pairIdx) = unitData.organoidLabel(targetRow);
        summarySourceNSpikes(pairIdx) = unitData.nSpikes(sourceRow);
        summaryTargetNSpikes(pairIdx) = unitData.nSpikes(targetRow);
        peakLatencyMs(pairIdx) = expectedLatenciesMs(peakIdx);
        peakCount(pairIdx) = thisPeak;
        baselineMeanCount(pairIdx) = thisBaseline;
        peakToBaselineRatio(pairIdx) = thisRatio;
        peakEfficacy(pairIdx) = thisEfficacy;
        nearZeroLagCount(pairIdx) = thisNearZero;
        sharedElectrodes(pairIdx) = thisShared;
        centroidDistanceUm(pairIdx) = thisDistance;
        nearZeroLagRatio(pairIdx) = thisNearZeroRatio;
        passesEfficacyCap(pairIdx) = passesEfficacy;
        passesDuplicateTest(pairIdx) = passesDuplicate;
        passesNearZeroLag(pairIdx) = passesNearZero;
        passesMinPeakCount(pairIdx) = passesCount;
        passesPeakToBaseline(pairIdx) = passesRatio;
        isConnection(pairIdx) = passes;

        rows = countRow + (1:nLatencyBins);
        sourceUnit(rows) = unitData.unitId(sourceRow);
        targetUnit(rows) = unitData.unitId(targetRow);
        sourceOrganoid(rows) = unitData.organoidLabel(sourceRow);
        targetOrganoid(rows) = unitData.organoidLabel(targetRow);
        latencyMs(rows) = expectedLatenciesMs;
        count(rows) = thisCounts;
        countRow = countRow + nLatencyBins;
    end
end

pairCounts = table(sourceUnit(1:countRow), targetUnit(1:countRow), ...
    sourceOrganoid(1:countRow), targetOrganoid(1:countRow), ...
    latencyMs(1:countRow), count(1:countRow), ...
    'VariableNames', {'sourceUnit', 'targetUnit', 'sourceOrganoid', ...
    'targetOrganoid', 'latencyMs', 'count'});

pairSummary = table(summarySource(1:pairIdx), summaryTarget(1:pairIdx), ...
    summarySourceOrganoid(1:pairIdx), summaryTargetOrganoid(1:pairIdx), ...
    summarySourceNSpikes(1:pairIdx), summaryTargetNSpikes(1:pairIdx), ...
    peakLatencyMs(1:pairIdx), peakCount(1:pairIdx), baselineMeanCount(1:pairIdx), ...
    peakToBaselineRatio(1:pairIdx), peakEfficacy(1:pairIdx), ...
    nearZeroLagCount(1:pairIdx), sharedElectrodes(1:pairIdx), ...
    centroidDistanceUm(1:pairIdx), nearZeroLagRatio(1:pairIdx), ...
    passesMinPeakCount(1:pairIdx), ...
    passesPeakToBaseline(1:pairIdx), passesEfficacyCap(1:pairIdx), ...
    passesDuplicateTest(1:pairIdx), passesNearZeroLag(1:pairIdx), ...
    isConnection(1:pairIdx), ...
    'VariableNames', {'sourceUnit', 'targetUnit', 'sourceOrganoid', ...
    'targetOrganoid', 'sourceNSpikes', 'targetNSpikes', 'peakLatencyMs', ...
    'peakCount', 'baselineMeanCount', 'peakToBaselineRatio', 'peakEfficacy', ...
    'nearZeroLagCount', 'sharedElectrodes', 'centroidDistanceUm', ...
    'nearZeroLagRatio', 'passesMinPeakCount', 'passesPeakToBaseline', ...
    'passesEfficacyCap', 'passesDuplicateTest', 'passesNearZeroLag', ...
    'isConnection'});
end
