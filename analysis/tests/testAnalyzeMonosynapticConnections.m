function tests = testAnalyzeMonosynapticConnections
% Synthetic regression tests for clustering, direction, and count thresholding.
tests = functiontests(localfunctions);
end

function testPositiveLagMeansTargetAfterSource(testCase)
sourceTimes = (1:20).';
targetTimes = sourceTimes + 0.005;
electrodeData = makeTwoUnitElectrodeTable(sourceTimes, targetTimes);

results = runAnalysis(electrodeData, 22, 5);

testCase.verifyEqual(height(results.unitData), 2);
testCase.verifyEqual(results.unitData.organoidLabel, ...
    ["Organoid 1"; "Organoid 2"]);
testCase.verifyEqual(height(results.connectionTable), 1);

connection = results.connectionTable(1, :);
testCase.verifyEqual(connection.sourceUnit, 1);
testCase.verifyEqual(connection.targetUnit, 2);
testCase.verifyEqual(connection.peakLatencyMs, 5, 'AbsTol', 1e-12);
testCase.verifyEqual(connection.peakCount, 20);
testCase.verifyEqual(connection.baselineMeanCount, 0);
testCase.verifyEqual(connection.peakToBaselineRatio, Inf);

summary = results.connectionSummaryTable;
testCase.verifyEqual(summary.relationship, ["Organoid 1 within"; ...
    "Organoid 2 within"; "Organoid 1 to Organoid 2"; ...
    "Organoid 2 to Organoid 1"]);
testCase.verifyEqual(summary.nConnections, [0; 0; 1; 0]);

reverse = results.pairSummaryTable( ...
    results.pairSummaryTable.sourceUnit == 2 & ...
    results.pairSummaryTable.targetUnit == 1, :);
testCase.verifyFalse(reverse.isConnection);
testCase.verifyEqual(reverse.peakCount, 0);
end

function testPeakToBaselineMultiplierIsLiteral(testCase)
baselineLagsMs = [(-11:0.5:-1), (5.5:0.5:11)].';
lagsSec = [baselineLagsMs / 1000; repmat(0.003, 5, 1)];
sourceTimes = (1:numel(lagsSec)).';
targetTimes = sourceTimes + lagsSec;
electrodeData = makeTwoUnitElectrodeTable(sourceTimes, targetTimes);

atFive = runAnalysis(electrodeData, 40, 5);
forwardAtFive = selectForwardPair(atFive);
testCase.verifyEqual(forwardAtFive.peakCount, 5);
testCase.verifyEqual(forwardAtFive.peakLatencyMs, 3, 'AbsTol', 1e-12);
testCase.verifyEqual(forwardAtFive.baselineMeanCount, 1, 'AbsTol', 1e-12);
testCase.verifyEqual(forwardAtFive.peakToBaselineRatio, 5, 'AbsTol', 1e-12);
testCase.verifyTrue(forwardAtFive.isConnection);

atSix = runAnalysis(electrodeData, 40, 6);
forwardAtSix = selectForwardPair(atSix);
testCase.verifyFalse(forwardAtSix.isConnection);
end

function testFiringFilterPrecedesDbscan(testCase)
sourceTimes = (1:20).';
targetTimes = sourceTimes + 0.005;
electrodeData = makeTwoUnitElectrodeTable(sourceTimes, targetTimes);
extraElectrode = table(99, 2000, 2000, {1}, 'VariableNames', ...
    electrodeData.Properties.VariableNames);
electrodeData = [electrodeData; extraElectrode];

results = runAnalysis(electrodeData, 22, 5);

testCase.verifyEqual(nnz(results.electrodeData.passesFiringFilter), 4);
testCase.verifyFalse(results.electrodeData.passesFiringFilter(end));
testCase.verifyEqual(height(results.unitData), 2);
end

function testUnassignedUnitsAreNotCorrelated(testCase)
sourceTimes = (1:20).';
targetTimes = sourceTimes + 0.005;
electrodeData = makeTwoUnitElectrodeTable(sourceTimes, targetTimes);
outsideUnit = table([31; 32], [2000; 2010], [2000; 2000], ...
    {sourceTimes; sourceTimes + 0.0002}, 'VariableNames', ...
    electrodeData.Properties.VariableNames);
electrodeData = [electrodeData; outsideUnit];

results = runAnalysis(electrodeData, 22, 5);

testCase.verifyEqual(height(results.unitData), 3);
testCase.verifyEqual(results.unitData.organoidLabel(3), "unassigned");
testCase.verifyEqual(height(results.pairSummaryTable), 2);
testCase.verifyFalse(any(results.pairSummaryTable.sourceUnit == 3));
testCase.verifyFalse(any(results.pairSummaryTable.targetUnit == 3));
end

function testBothUnitsMustHaveStrictlyMoreThanFiftySpikes(testCase)
fiftySource = (1:50).';
fiftyTarget = fiftySource + 0.005;
atThreshold = makeTwoUnitElectrodeTable(fiftySource, fiftyTarget);
atThresholdResults = runAnalysis(atThreshold, 52, 5, 50);

testCase.verifyFalse(any(atThresholdResults.unitData.passesPairSpikeFilter));
testCase.verifyEqual(height(atThresholdResults.pairSummaryTable), 0);

fiftyOneSource = (1:51).';
fiftyOneTarget = fiftyOneSource + 0.005;
aboveThreshold = makeTwoUnitElectrodeTable(fiftyOneSource, fiftyOneTarget);
aboveThresholdResults = runAnalysis(aboveThreshold, 53, 5, 50);

testCase.verifyTrue(all(aboveThresholdResults.unitData.passesPairSpikeFilter));
testCase.verifyEqual(height(aboveThresholdResults.pairSummaryTable), 2);
forward = selectForwardPair(aboveThresholdResults);
testCase.verifyEqual(forward.sourceNSpikes, 51);
testCase.verifyEqual(forward.targetNSpikes, 51);
testCase.verifyTrue(forward.isConnection);
end

function testMinimumPeakCountRejectsSingleCoincidence(testCase)
sourceTimes = (1:51).';
targetTimes = sourceTimes + 0.5;
targetTimes(1) = sourceTimes(1) + 0.005;
electrodeData = makeTwoUnitElectrodeTable(sourceTimes, targetTimes);

results = runAnalysis(electrodeData, 53, 5, 50);
forward = selectForwardPair(results);

testCase.verifyEqual(forward.sourceNSpikes, 51);
testCase.verifyEqual(forward.targetNSpikes, 51);
testCase.verifyEqual(forward.peakCount, 1);
testCase.verifyEqual(forward.baselineMeanCount, 0);
testCase.verifyTrue(forward.passesPeakToBaseline);
testCase.verifyFalse(forward.passesMinPeakCount);
testCase.verifyFalse(forward.isConnection);
end

function testTenMillisecondPeakIsIncludedInCausalWindow(testCase)
sourceTimes = (1:51).';
targetTimes = sourceTimes + 0.010;
electrodeData = makeTwoUnitElectrodeTable(sourceTimes, targetTimes);

results = runAnalysis(electrodeData, 53, 5, 50);
forward = selectForwardPair(results);

testCase.verifyEqual(forward.peakCount, 51);
testCase.verifyEqual(forward.peakLatencyMs, 10, 'AbsTol', 1e-12);
testCase.verifyTrue(forward.passesMinPeakCount);
testCase.verifyTrue(forward.isConnection);
end

function testFlatCorrelogramIsRejectedDespitePeakCount(testCase)
allLagsMs = [(-11:0.5:-1), (1:0.5:11)].';
lagsSec = repelem(allLagsMs, 10) / 1000;
sourceTimes = (1:numel(lagsSec)).';
targetTimes = sourceTimes + lagsSec;
electrodeData = makeTwoUnitElectrodeTable(sourceTimes, targetTimes);

results = runAnalysis(electrodeData, numel(lagsSec) + 2, 5);
forward = selectForwardPair(results);

testCase.verifyEqual(forward.peakCount, 10);
testCase.verifyEqual(forward.baselineMeanCount, 10, 'AbsTol', 1e-12);
testCase.verifyTrue(forward.passesMinPeakCount);
testCase.verifyFalse(forward.passesPeakToBaseline);
testCase.verifyFalse(forward.isConnection);
end

function results = runAnalysis(electrodeData, durationSec, multiplier, minUnitSpikes)
if nargin < 4
    minUnitSpikes = 0;
end
results = analyze_monosynaptic_connections(electrodeData, ...
    'RecordingDurationSec', durationSec, ...
    'MinElectrodeFiringRateHz', 0.1, ...
    'MaxElectrodeFiringRateHz', Inf, ...
    'DbscanRadiusUm', 50, ...
    'DbscanMinElectrodes', 2, ...
    'MergeWindowMs', 1, ...
    'MinUnitSpikesForPair', minUnitSpikes, ...
    'CorrelationBinMs', 0.5, ...
    'LatencyRangeMs', [1 11], ...
    'PeakLatencyRangeMs', [1 11], ...
    'MinPeakCount', 5, ...
    'PeakToBaselineMultiplier', multiplier, ...
    'OrganoidLabels', {'Organoid 1', 'Organoid 2'}, ...
    'OrganoidRoiPositions', [-50 -50 200 100; 900 -50 200 100], ...
    'InteractiveRoiSelection', false);
end

function pair = selectForwardPair(results)
pair = results.pairSummaryTable( ...
    results.pairSummaryTable.sourceUnit == 1 & ...
    results.pairSummaryTable.targetUnit == 2, :);
end

function electrodeData = makeTwoUnitElectrodeTable(sourceTimes, targetTimes)
electrode = [11; 12; 21; 22];
x = [0; 10; 1000; 1010];
y = [0; 0; 0; 0];
spikeTimesSec = {sourceTimes; sourceTimes + 0.0002; ...
    targetTimes; targetTimes + 0.0002};
electrodeData = table(electrode, x, y, spikeTimesSec);
end
