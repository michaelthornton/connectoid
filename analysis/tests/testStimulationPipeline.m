function tests = testStimulationPipeline
% The stimulation analysis end to end on a synthetic connectoid: splitting the
% organoids, both metrics, the summary row, the running table and the figures.
% Everything downstream of reading the H5.
tests = functiontests(localfunctions);
end


function setupOnce(testCase)
here = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(here, 'shared'), fullfile(here, 'stim'));
testCase.TestData.work = tempname;
mkdir(testCase.TestData.work);
end


function teardownOnce(testCase)
rmdir(testCase.TestData.work, 's');
close all;
end


function testRespondingDirectionIsCalledConnected(testCase)
rois = roiTable();
withResponse = synthetic("DG", 0.45, 1);
withoutResponse = synthetic("CA3", 0, 2);
recs = [withResponse, withoutResponse];

clear metrics
for p = 1:2
    undriven = stim_organoids(recs(p), rois);
    metrics(p) = stim_population_metrics(recs(p), undriven, ...
        'Iterations', 2000, 'Seed', 1); %#ok<SAGROW>
end

% Recording 1 stimulated DG, and CA3 was given an evoked response.
testCase.verifyEqual(metrics(1).respondingOrganoid, "CA3");
testCase.verifyEqual(metrics(2).respondingOrganoid, "DG");
testCase.verifyTrue(metrics(1).ciSeparated);
testCase.verifyFalse(metrics(2).ciSeparated);
testCase.verifyGreaterThan(metrics(1).fracAbove, 20);
% fracAbove has a no-effect value of 5, so the quiet direction sits near it.
testCase.verifyLessThan(metrics(2).fracAbove, 20);
end


function testSummaryRowAndRunningTable(testCase)
rois = roiTable();
recs = [synthetic("DG", 0.45, 1), synthetic("CA3", 0, 2)];
clear metrics
for p = 1:2
    metrics(p) = stim_population_metrics(recs(p), stim_organoids(recs(p), rois), ...
        'Iterations', 2000, 'Seed', 1); %#ok<SAGROW>
end

row = stim_summary_row("synthetic01", recs, rois, metrics);
testCase.verifyEqual(height(row), 1);
testCase.verifyEqual(row.stimA, "DG");
testCase.verifyEqual(row.stimA_respondingOrganoid, "CA3");
testCase.verifyEqual(row.stimA_connected, 1);
testCase.verifyEqual(row.stimB_connected, 0);
testCase.verifyEqual(row.outcome, "DG->CA3 only");

file = fullfile(testCase.TestData.work, "table.csv");
append_summary_row(file, row);
append_summary_row(file, row);
% Analysing the same connectoid twice leaves one row.
T = readtable(file, 'TextType', 'string');
testCase.verifyEqual(height(T), 1);

second = row;
second.dataset = "synthetic02";
append_summary_row(file, second);
testCase.verifyEqual(height(readtable(file)), 2);
end


function testBothFiguresDraw(testCase)
rois = roiTable();
recs = [synthetic("DG", 0.45, 1), synthetic("CA3", 0, 2)];
clear metrics
for p = 1:2
    metrics(p) = stim_population_metrics(recs(p), stim_organoids(recs(p), rois), ...
        'Iterations', 500, 'Seed', 1); %#ok<SAGROW>
end
f1 = stim_plot_dataset(recs, rois, metrics, "synthetic01");
testCase.verifyTrue(isgraphics(f1, 'figure'));
f2 = stim_plot_response_evidence(metrics, ["DG -> CA3", "CA3 -> DG"], "synthetic01");
testCase.verifyTrue(isgraphics(f2, 'figure'));
close([f1 f2]);
end


function R = roiTable()
% Two boxes far enough apart that no electrode falls in both.
R = table(["DG"; "CA3"], [0; 2000], [0; 0], [600; 600], [600; 600], ...
    'VariableNames', {'type', 'x', 'y', 'width', 'height'});
end


function S = synthetic(stimulated, responseProb, seed)
% One stimulation recording: Poisson baseline on both organoids, trains once a
% second, and an evoked burst in the undriven organoid on responseProb of
% the trials.
rng(seed, 'twister');
dg = (0:9)';            % electrode ids, placed inside the DG box below
ca3 = (100:109)';
S.recording = "synthetic_" + stimulated;
S.file = "synthetic/" + stimulated + "/data.raw.h5";
S.electrodes = [dg; ca3];
S.x = [linspace(50, 550, 10)'; linspace(2050, 2550, 10)'];
S.y = [linspace(50, 550, 10)'; linspace(50, 550, 10)'];
S.stimElectrodes = dg;
if stimulated == "CA3", S.stimElectrodes = ca3; end

baselineSec = 60;
nTrials = 60;
S.onsets = baselineSec + (0:nTrials - 1)';
S.duration = S.onsets(end) + 1;

rateHz = 4;
t = []; el = [];
for e = [dg; ca3]'
    n = poissrnd(rateHz * S.duration);
    t = [t; S.duration * rand(n, 1)]; %#ok<AGROW>
    el = [el; repmat(e, n, 1)]; %#ok<AGROW>
end

% Extra spikes 30 to 60 ms after the onset, in the undriven organoid.
responders = ca3;
if stimulated == "CA3", responders = dg; end
for k = 1:nTrials
    if rand > responseProb, continue; end
    for e = responders'
        n = poissrnd(3);
        t = [t; S.onsets(k) + 0.030 + 0.030 * rand(n, 1)]; %#ok<AGROW>
        el = [el; repmat(e, n, 1)]; %#ok<AGROW>
    end
end
[S.t, order] = sort(t);
S.el = el(order);
S.amp = ones(size(S.t));
end
