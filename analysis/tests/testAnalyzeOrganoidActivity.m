function tests = testAnalyzeOrganoidActivity
% Synthetic regression tests for unit firing rates and MFR peak selection.
tests = functiontests(localfunctions);
end

function testUnitMeansAndDetectedPeaks(testCase)
durationSec = 10;
peakStarts = [1; 4; 7];
organoid1TimesA = sort([peakStarts; peakStarts + 0.02; 0.3; 2.3; 5.3; 8.3]);
organoid1TimesB = sort([peakStarts + 0.01; peakStarts + 0.03; 0.7; 2.7; 5.7; 8.7]);
organoid2TimesA = (0.25:1:9.25).';
organoid2TimesB = (0.75:1:9.75).';

unitData = table((1:4).', ...
    {organoid1TimesA; organoid1TimesB; organoid2TimesA; organoid2TimesB}, ...
    ["Organoid 1"; "Organoid 1"; "Organoid 2"; "Organoid 2"], ...
    'VariableNames', {'unitId', 'spikeTimesSec', 'organoidLabel'});

activity = analyze_organoid_activity(unitData, ...
    'RecordingDurationSec', durationSec, ...
    'RateBinSec', 0.01, ...
    'RateGaussianSigmaSec', 0.03, ...
    'PeakMinimumDistanceSec', 2, ...
    'PeakMinimumProminenceHz', 1);

summary = activity.summaryTable;
testCase.verifyEqual(summary.nUnits, [2; 2]);
testCase.verifyEqual(summary.meanUnitFiringRateHz, [10/10; 10/10], ...
    'AbsTol', 1e-12);
testCase.verifyGreaterThanOrEqual(summary.peakCount(1), 3);
testCase.verifyGreaterThan(summary.meanPeakFiringRateHz(1), 1);
testCase.verifyGreaterThan(summary.peakToMeanRateRatio(1), 1);
testCase.verifyEqual(height(activity.peakTable), sum(summary.peakCount));
testCase.verifySize(activity.unitRateHz, [numel(activity.traceTimeSec), 4]);
testCase.verifySize(activity.timingRateHz, [numel(activity.traceTimeSec), 2]);
testCase.verifyEqual(numel(activity.timingPeakTimesSec{1}), summary.peakCount(1));
testCase.verifyLessThanOrEqual(max(abs(activity.timingPeakTimesSec{1} - ...
    activity.peakTimesSec{1})), 0.5 + 1e-12);
end

function testMinimumDistanceSeparatesDetectedPeaks(testCase)
durationSec = 8;
spikes = sort([1; 1.01; 1.5; 1.51; 5; 5.01; 0.2; 2.8; 4; 7]);
unitData = table(1, {spikes}, "Organoid 1", ...
    'VariableNames', {'unitId', 'spikeTimesSec', 'organoidLabel'});

activity = analyze_organoid_activity(unitData, ...
    'RecordingDurationSec', durationSec, ...
    'RateBinSec', 0.01, ...
    'RateGaussianSigmaSec', 0.02, ...
    'PeakMinimumDistanceSec', 2, ...
    'PeakMinimumProminenceHz', 0.1);

times = activity.peakTimesSec{1};
testCase.verifyGreaterThanOrEqual(numel(times), 2);
testCase.verifyGreaterThanOrEqual(min(diff(times)), 2 - 1e-12);
testCase.verifyEqual(activity.parameters.PeakMinimumDistanceSec, 2);
end
