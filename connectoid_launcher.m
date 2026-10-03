%% DG-CA3 connectoid launcher
% One spontaneous recording of a two-organoid connectoid. Clusters spatial
% units, then monosynaptic connectivity, unit activity and inter-organoid
% synchrony. Run the sections in order.
%
% Stimulation recordings: stim_launcher.m.

clear; close all; clc;
set(groot,'defaultFigureWindowStyle','docked');
set(groot,'defaultFigureRenderer','painters', ...
    'defaultAxesFontName','Helvetica','defaultTextFontName','Helvetica');

%% Settings

config = struct();
config.DataFile = ''; %path to network scan .h5 file or folder
config.WellID = 1;
config.RecordingNumber = 1;

% Electrode firing-rate filter and spatial unit clustering.
config.MinElectrodeFiringRateHz = 0.01;
config.MaxElectrodeFiringRateHz = Inf;
config.DbscanRadiusUm = 60;
config.DbscanMinElectrodes = 2;
config.MergeWindowMs = 1.5;

config.MinUnitSpikesForPair = 50;    % both units of a pair need more than this

% Connection rule. Each of these changes which pairs are accepted.
config.CorrelationBinMs = 0.5;
config.LatencyRangeMs = [1 11];      % counts retained over -11 to +11 ms
config.PeakLatencyRangeMs = [1 11];  % causal peak searched here only
config.MinPeakCount = 10;
config.PeakToBaselineMultiplier = 5;
config.MaxPeakEfficacy = 0.5;        % source spikes followed by a target spike
config.DuplicateRadiusUm = 60;       % units closer than this are one neuron
config.MaxNearZeroLagRatio = 2;      % |lag| < 1 ms per bin over the peak bin

% Unit activity and peaks in the organoid mean firing rate.
config.FiringRateBinSec = 0.01;
config.FiringRateGaussianSigmaSec = 0.5;
config.SynchronyTimingGaussianSigmaSec = 0.05;
config.SynchronyTimingRefinementWindowSec = 0.5;
config.MfrPeakMinimumDistanceSec = 2;
config.MfrPeakMinimumProminenceHz = 1;

% Inter-organoid timing, against a circular-shift null.
config.SynchronyDelayBinSec = 0.05;
config.SynchronyMaxLagSec = 1.0;
config.SynchronyNumCircularShifts = 500;
config.SynchronyRandomSeed = 1;

% Organoid assignment. The ROI step is where it happens: each box is drawn and
% labelled with its type, and everything downstream carries that name. Draw them
% in this order. Leave positions empty to draw the ROIs, then paste the printed
% positions back here for a non-interactive rerun.
config.OrganoidTypes = {'DG','CA3'};
config.OrganoidRoiPositions = [];      % N-by-4 rows: [x y width height]
config.UseSquareRois = true;

config.MaxCorrelogramsToPlot = 12;
config.OutputDirectory = "";           % empty writes beside the data file
config.Label = "";                     % names the row; empty uses the recording id
config.TableFile = "";                 % empty writes connectoid_spontaneous.csv

%% Paths

% Running a single section executes a temporary copy of this file, so mfilename
% alone is not enough to locate the project.
candidates = {fileparts(mfilename('fullpath')), pwd, fileparts(which('connectoid_launcher'))};
root = '';
for k = 1:numel(candidates)
    if ~isempty(candidates{k}) && isfolder(fullfile(candidates{k}, 'analysis'))
        root = candidates{k};
        break
    end
end
assert(~isempty(root), 'connectoid:ProjectRoot', ...
    ['Could not find the DG_CA3_Analysis folder. Make it MATLAB''s Current ' ...
     'Folder, or add it to the path, then run this section again.' ...
     '\nCurrent folder: %s'], pwd);
addpath(fullfile(root, 'analysis', 'shared'), ...
        fullfile(root, 'analysis', 'spontaneous'));
add_maxwell_path(root);

if strlength(config.DataFile)==0
    [f,p] = uigetfile({'*.h5;*.raw.h5;*.data.raw.h5','MaxWell H5 recordings'; ...
        '*.*','All files'},'Select a MaxWell network recording');
    if isequal(f,0), error('connectoid:Cancelled','No recording selected.'); end
    config.DataFile = string(fullfile(p,f));
end
config.DataFile = resolve_recording(config.DataFile);

%% Monosynaptic connectivity

fprintf('\n=== Monosynaptic connectivity analysis ===\n');
if isempty(config.OrganoidRoiPositions)
    fprintf('Draw one box per organoid in the window that opens, then Confirm.\n');
end
results = analyze_monosynaptic_connections(config.DataFile, ...
    'WellID',config.WellID,'RecordingNumber',config.RecordingNumber, ...
    'MinElectrodeFiringRateHz',config.MinElectrodeFiringRateHz, ...
    'MaxElectrodeFiringRateHz',config.MaxElectrodeFiringRateHz, ...
    'DbscanRadiusUm',config.DbscanRadiusUm, ...
    'DbscanMinElectrodes',config.DbscanMinElectrodes, ...
    'MergeWindowMs',config.MergeWindowMs, ...
    'MinUnitSpikesForPair',config.MinUnitSpikesForPair, ...
    'CorrelationBinMs',config.CorrelationBinMs, ...
    'LatencyRangeMs',config.LatencyRangeMs, ...
    'PeakLatencyRangeMs',config.PeakLatencyRangeMs, ...
    'MinPeakCount',config.MinPeakCount, ...
    'PeakToBaselineMultiplier',config.PeakToBaselineMultiplier, ...
    'MaxPeakEfficacy',config.MaxPeakEfficacy, ...
    'DuplicateRadiusUm',config.DuplicateRadiusUm, ...
    'MaxNearZeroLagRatio',config.MaxNearZeroLagRatio, ...
    'OrganoidLabels',config.OrganoidTypes, ...
    'OrganoidRoiPositions',config.OrganoidRoiPositions, ...
    'InteractiveRoiSelection',isempty(config.OrganoidRoiPositions), ...
    'UseSquareRois',config.UseSquareRois);

P = results.parameters;
fprintf('Electrodes passing firing filter: %d / %d\n', ...
    nnz(results.electrodeData.passesFiringFilter),height(results.electrodeData));
fprintf('DBSCAN units: %d (assigned to an organoid: %d)\n', ...
    height(results.unitData),nnz(results.unitData.organoidLabel~="unassigned"));
fprintf('Assigned units with > %d spikes: %d\n',P.MinUnitSpikesForPair, ...
    nnz(results.unitData.organoidLabel~="unassigned" & results.unitData.passesPairSpikeFilter));
fprintf('Directed pairs evaluated: %d\n',height(results.pairSummaryTable));
fprintf(['Connection rule: peak at +%.3g to +%.3g ms, peak >= %d counts, ' ...
    'and peak >= %.3g x mean off-peak baseline\n'], ...
    P.PeakLatencyRangeMs(1),P.PeakLatencyRangeMs(2), ...
    P.MinPeakCount,P.PeakToBaselineMultiplier);
fprintf('Connections passing both rules: %d\n',height(results.connectionTable));
disp(results.connectionTable);
disp(results.connectionSummaryTable);

%% Unit activity

fprintf('\n=== Organoid unit activity analysis ===\n');
results.activity = analyze_organoid_activity(results.unitData, ...
    'RecordingDurationSec',P.RecordingDurationSec, ...
    'OrganoidLabels',config.OrganoidTypes, ...
    'RateBinSec',config.FiringRateBinSec, ...
    'RateGaussianSigmaSec',config.FiringRateGaussianSigmaSec, ...
    'TimingGaussianSigmaSec',config.SynchronyTimingGaussianSigmaSec, ...
    'TimingRefinementWindowSec',config.SynchronyTimingRefinementWindowSec, ...
    'PeakMinimumDistanceSec',config.MfrPeakMinimumDistanceSec, ...
    'PeakMinimumProminenceHz',config.MfrPeakMinimumProminenceHz);
fprintf('MFR peak rule: >= %.3g Hz local prominence, >= %.3g s apart\n', ...
    config.MfrPeakMinimumProminenceHz,config.MfrPeakMinimumDistanceSec);
disp(results.activity.summaryTable);

%% Inter-organoid synchrony

fprintf('\n=== Inter-organoid synchrony analysis ===\n');
results.synchrony = analyze_organoid_synchrony(results.activity, ...
    'DelayBinSec',config.SynchronyDelayBinSec, ...
    'MaxLagSec',config.SynchronyMaxLagSec, ...
    'NumCircularShifts',config.SynchronyNumCircularShifts, ...
    'RandomSeed',config.SynchronyRandomSeed);
disp(results.synchrony.summaryTable);
disp(results.synchrony.directionalSummaryTable);

%% Figures and export

out = resolveOutputDirectory(config);
if ~isfolder(out)
    [ok,msg] = mkdir(out);
    assert(ok,'connectoid:OutputDirectory','Could not create %s\n%s',out,msg);
end

figures = plot_connectoid_overview(results,config.MaxCorrelogramsToPlot);
save_figure_set(figures,out,["01_electrode_firing_filter";"02_dbscan_units_and_organoids"; ...
    "03_directed_connections";"04_connection_counts"; ...
    "05_organoid_activity";"06_interorganoid_synchrony"]);

%% Summary row

if strlength(config.Label) == 0
    config.Label = recording_label(config.DataFile);
end
row = connectoid_summary_row(config.Label, results, string(config.OrganoidTypes));
disp(row);
tableFile = config.TableFile;
if strlength(tableFile) == 0
    tableFile = fullfile(fileparts(fileparts(out)), "connectoid_spontaneous.csv");
end
append_summary_row(tableFile, row);

%% Save

config.OrganoidRoiPositions = results.organoidRois.position;
results.launcherConfig = config;
save(fullfile(out,'experimental_analysis_results.mat'),'results','config','-v7.3');
export_connectoid_results(results,out);

snapshot_source(root,out,'connectoid_launcher.m', ...
    ["analyze_monosynaptic_connections","analyze_organoid_activity", ...
     "analyze_organoid_synchrony"]);

fprintf('\nSaved reproducible analysis to:\n%s\n',out);
fprintf('For a non-interactive rerun, use these ROI positions:\n');
disp(results.organoidRois);

%% Local functions

function out = resolveOutputDirectory(config)
if strlength(config.OutputDirectory)>0
    out = char(config.OutputDirectory); return;
end
[dataPath,dataName] = fileparts(char(config.DataFile));
out = fullfile(dataPath,[regexprep(dataName,'[^A-Za-z0-9_-]','_') '_analysis']);
end

function snapshot_source(root,out,launcherName,functionNames)
% Copy the code that produced this result next to it.
source = fullfile(root,launcherName);
if isfile(source)
    copyfile(source,fullfile(out,[extractBefore(launcherName,'.m') '_snapshot.m']));
end
for name = functionNames
    found = which(char(name));
    if ~isempty(found)
        copyfile(found,fullfile(out,char(name)+"_snapshot.m"));
    end
end
end
