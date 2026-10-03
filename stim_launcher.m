%% Connectoid stimulation launcher
% Two MaxWell stimulation recordings of one pair of organoids, one with each
% organoid stimulated. Run the sections in order.

clear; close all; clc;
set(groot, 'defaultFigureWindowStyle', 'docked');
set(groot, 'defaultFigureRenderer', 'painters', ...
    'defaultAxesFontName', 'Helvetica', 'defaultTextFontName', 'Helvetica');

%% Settings

config = struct();

% Each recording is either a data.raw.h5 or the folder holding it. Leave both
% empty to choose them in a dialog. OrganoidTypes(1) is the organoid that was
% stimulated in StimulateFirst.
config.StimulateFirst = ''; %path to first stim file or folder
config.StimulateSecond = ''; %path to second stim file or folder
config.OrganoidTypes = ["DG" "CA3"];
config.WellID = 1;
config.Label = "";                  % names the row; empty uses the recording id

% Organoid boundaries in micrometres. Leave both empty to draw them, then reuse
% them on later runs through either setting.
config.OrganoidRoiFile = ''; %path to ROI file if already ran connectoid_launcher
config.OrganoidRoiPositions = [];   % 2-by-4: [x y width height] per organoid
config.UseSquareRois = true;

% Response window and the blank before it. WindowSec(1) must clear the last
% pulse of the train.
config.WindowSec = [0.030 0.100];
config.PreBlankSec = 0.005;
config.Iterations = 20000;
config.GridStepSec = 0.001;
config.Seed = 1;

% Figures. Gamma and SmoothBins affect display only.
config.SpanSec = [-0.100 0.400];
config.BinSec = 0.020;              % heatmap bins
config.CurveBinSec = 0.010;         % rate panel bins
config.ScalePct = 99;
config.Gamma = 0.6;
config.SmoothBins = 0.45;

config.OutputDirectory = "";        % empty writes beside the recordings
config.TableFile = "";              % empty writes beside the output directory

%% Paths

% Running a single section executes a temporary copy of this file, so mfilename
% alone is not enough to locate the project.
candidates = {fileparts(mfilename('fullpath')), pwd, fileparts(which('stim_launcher'))};
root = '';
for k = 1:numel(candidates)
    if ~isempty(candidates{k}) && isfolder(fullfile(candidates{k}, 'analysis'))
        root = candidates{k};
        break
    end
end
assert(~isempty(root), 'stim:ProjectRoot', ...
    ['Could not find the DG_CA3_Analysis folder. Make it MATLAB''s Current ' ...
     'Folder, or add it to the path, then run this section again.' ...
     '\nCurrent folder: %s'], pwd);
addpath(fullfile(root, 'analysis', 'shared'), fullfile(root, 'analysis', 'stim'));
add_maxwell_path(root);

%% Select the recordings

for k = ["StimulateFirst" "StimulateSecond"]
    field = char(k);
    if strlength(config.(field)) > 0, continue; end
    stimulated = config.OrganoidTypes(1 + (k == "StimulateSecond"));
    [f, folder] = uigetfile({'*.h5', 'MaxWell H5 recordings'; '*.*', 'All files'}, ...
        sprintf('Select the recording where %s was stimulated', stimulated));
    if isequal(f, 0), error('stim:Cancelled', 'No recording selected.'); end
    config.(field) = string(fullfile(folder, f));
end
config.StimulateFirst = resolve_recording(config.StimulateFirst);
config.StimulateSecond = resolve_recording(config.StimulateSecond);
files = [config.StimulateFirst config.StimulateSecond];
assert(files(1) ~= files(2), 'stim:SameFile', ...
    'Both directions point at the same recording.');
if strlength(config.Label) == 0
    config.Label = recording_label(files(1));
end
fprintf('Connectoid %s\n  %s stimulated: %s\n  %s stimulated: %s\n', ...
    config.Label, config.OrganoidTypes(1), files(1), config.OrganoidTypes(2), files(2));

%% Load

clear recs
for p = 1:2
    recs(p) = stim_read_recording(files(p), 'WellID', config.WellID); %#ok<SAGROW>
    fprintf('%s: %d routed electrodes, %d trains, %.0f s\n', ...
        recs(p).recording, numel(recs(p).electrodes), numel(recs(p).onsets), ...
        recs(p).duration);
end

%% Organoid boundaries

if strlength(config.OrganoidRoiFile) > 0
    rois = stim_read_rois(config.OrganoidRoiFile);
else
    if isempty(config.OrganoidRoiPositions)
        fprintf('Draw one box per organoid in the window that opens, then Confirm.\n');
        config.OrganoidRoiPositions = select_organoid_rois( ...
            electrodeTable(recs), config.OrganoidTypes(:), config.UseSquareRois);
    end
    assert(isequal(size(config.OrganoidRoiPositions), [2 4]), 'stim:RoiSize', ...
        'OrganoidRoiPositions needs one [x y width height] row per organoid.');
    rois = table(config.OrganoidTypes(:), config.OrganoidRoiPositions(:, 1), ...
        config.OrganoidRoiPositions(:, 2), config.OrganoidRoiPositions(:, 3), ...
        config.OrganoidRoiPositions(:, 4), ...
        'VariableNames', {'type', 'x', 'y', 'width', 'height'});
end
disp(rois);

% Which organoid was stimulated is read from the assay, not from the boxes.
driven = strings(1, 2);
for p = 1:2
    [~, d] = stim_organoids(recs(p), rois);
    driven(p) = d.type;
end
assert(numel(unique(driven)) == 2, 'stim:SameDriven', ...
    ['Both recordings stimulated %s. Either the two files are the same ' ...
     'direction, or the boxes are in the wrong order.'], driven(1));

%% Response metrics

clear metrics
for p = 1:2
    [undriven, ~] = stim_organoids(recs(p), rois);
    metrics(p) = stim_population_metrics(recs(p), undriven, ...
        'WindowSec', config.WindowSec, 'Iterations', config.Iterations, ...
        'GridStepSec', config.GridStepSec, 'Seed', config.Seed); %#ok<SAGROW>
    fprintf('\n%s -> %s\n', driven(p), undriven.type);
    fprintf('  window mean %+.1f%% of baseline, 95%% CI +/- %.1f, p %.4f%s\n', ...
        metrics(p).deltaPct, metrics(p).deltaPctCI, metrics(p).pDelta, ...
        connectedNote(metrics(p).ciSeparated));
    fprintf('  %.1f%% of trials above the baseline 95th percentile, p %.4f\n', ...
        metrics(p).fracAbove, metrics(p).pFrac);
end

%% Summary row

row = stim_summary_row(config.Label, recs, rois, metrics, ...
    'WindowSec', config.WindowSec);
disp(row);

out = resolveOutputDirectory(config, files(1));
if ~isfolder(out)
    [ok, msg] = mkdir(out);
    assert(ok, 'stim:OutputDirectory', 'Could not create %s\n%s', out, msg);
end
tableFile = config.TableFile;
if strlength(tableFile) == 0
    tableFile = fullfile(fileparts(out), "connectoid_stimulation.csv");
end
groupTable = append_summary_row(tableFile, row);

%% Figures

fig1 = stim_plot_dataset(recs, rois, metrics, config.Label, ...
    'SpanSec', config.SpanSec, 'BinSec', config.BinSec, ...
    'CurveBinSec', config.CurveBinSec, 'WindowSec', config.WindowSec, ...
    'PreBlankSec', config.PreBlankSec, 'ScalePct', config.ScalePct, ...
    'Gamma', config.Gamma, 'SmoothBins', config.SmoothBins);
stim_save_figure(fig1, fullfile(out, "01_trials_and_rate_" + config.Label));

fig2 = stim_plot_response_evidence(metrics, ...
    driven + " -> " + [metrics.respondingOrganoid], config.Label, ...
    'WindowSec', config.WindowSec);
stim_save_figure(fig2, fullfile(out, "02_response_evidence_" + config.Label));

%% Save

config.OrganoidRoiPositions = [rois.x rois.y rois.width rois.height];
writetable(rois, fullfile(out, "organoid_rois.csv"));
save(fullfile(out, 'stimulation_analysis_results.mat'), ...
    'config', 'rois', 'metrics', 'row', '-v7.3');
copyfile(fullfile(root, 'stim_launcher.m'), fullfile(out, 'stim_launcher_snapshot.m'));

fprintf('\nSaved to:\n%s\n', out);
fprintf('Group table: %s (%d connectoids)\n', tableFile, height(groupTable));
fprintf('To rerun without drawing the boxes again, set:\n');
fprintf('  config.OrganoidRoiFile = "%s";\n', fullfile(out, "organoid_rois.csv"));

%% Local functions

function s = connectedNote(tf)
if tf
    s = '   CONNECTED';
else
    s = '';
end
end

function out = resolveOutputDirectory(config, firstFile)
if strlength(config.OutputDirectory) > 0
    out = string(config.OutputDirectory); return;
end
% Beside the recording folder, not inside it.
folder = fileparts(char(firstFile));
out = fullfile(fileparts(folder), char(config.Label) + "_stim_analysis");
end

function T = electrodeTable(recs)
% Routed electrodes and spike counts pooled over both recordings.
x = []; y = []; el = []; n = [];
for p = 1:numel(recs)
    S = recs(p);
    [tf, loc] = ismember(S.el, S.electrodes);
    counts = accumarray(loc(tf), 1, [numel(S.electrodes) 1]);
    el = [el; S.electrodes]; x = [x; S.x]; y = [y; S.y]; n = [n; counts]; %#ok<AGROW>
end
[uel, ~, g] = unique(el);
[~, first] = ismember(uel, el);
T = table(x(first), y(first), accumarray(g, n), ...
    'VariableNames', {'x', 'y', 'nSpikes'});
end
