function tests = testAnalyzeOrganoidSynchrony
% Synthetic regression tests for all-pairs inter-organoid event timing.
tests = functiontests(localfunctions);
end

function testPostTriggerPairsAndCircularNullAreReproducible(testCase)
sourceTimes = [1; 3; 5; 7; 9];
targetTimes = sourceTimes + 0.5;
activity = struct();
activity.organoidLabels = ["Organoid 1"; "Organoid 2"];
activity.parameters = struct('RecordingDurationSec', 10);
activity.peakTimesSec = {sourceTimes; targetTimes};

synchrony = analyze_organoid_synchrony(activity, ...
    'DelayBinSec', 0.25, 'MaxLagSec', 1, ...
    'NumCircularShifts', 20, 'RandomSeed', 11);
repeat = analyze_organoid_synchrony(activity, ...
    'DelayBinSec', 0.25, 'MaxLagSec', 1, ...
    'NumCircularShifts', 20, 'RandomSeed', 11);

pairs = synchrony.postTriggerPairTable;
organoidOnePairs = pairs(pairs.sourceOrganoid == "Organoid 1", :);
testCase.verifyEqual(height(organoidOnePairs), 5);
testCase.verifyEqual(organoidOnePairs.targetDelaySec, repmat(0.5, 5, 1), ...
    'AbsTol', 1e-12);
testCase.verifyEqual(synchrony.directionalSummaryTable.modalPairCount(1), 5);
testCase.verifyEqual(synchrony.directionalSummaryTable.modalPairCount(2), 0);
testCase.verifyEqual(synchrony.nullGlobalMaximumPairCount, ...
    repeat.nullGlobalMaximumPairCount, 'AbsTol', 1e-12);
end

function testEmptyPeakTrainsReturnZeroPairCounts(testCase)
activity = struct();
activity.organoidLabels = ["Organoid 1"; "Organoid 2"];
activity.parameters = struct('RecordingDurationSec', 1);
activity.peakTimesSec = {zeros(0, 1); zeros(0, 1)};

synchrony = analyze_organoid_synchrony(activity, ...
    'DelayBinSec', 0.1, 'MaxLagSec', 0.5, ...
    'NumCircularShifts', 5, 'RandomSeed', 1);
testCase.verifyEqual(height(synchrony.postTriggerPairTable), 0);
testCase.verifyEqual(synchrony.directionalSummaryTable.modalPairCount, [0; 0]);
testCase.verifyEqual(synchrony.nullGlobalMaximumPairCount, zeros(5, 1));
end
