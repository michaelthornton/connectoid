function row = connectoid_summary_row(dataset, results, types)
% One spontaneous recording as one row, for the running table.
%
% types names what each ROI is, in the order of the analysis labels.
arguments
    dataset (1,1) string
    results struct
    types (1,2) string
end
labels = string(results.parameters.OrganoidLabels(:))';
assert(numel(labels) == 2, 'connectoid:OrganoidCount', ...
    'Expected two organoids, got %d.', numel(labels));
act = results.activity.summaryTable;
con = results.connectionSummaryTable;
S = results.synchrony.summaryTable;

row = table(dataset, string(results.provenance.Source.DataFile), ...
    results.parameters.RecordingDurationSec, ...
    nnz(results.electrodeData.passesFiringFilter), height(results.unitData), ...
    'VariableNames', {'dataset', 'file', 'durationSec', ...
                      'nElectrodesPassingFilter', 'nUnits'});
row = [row, side("A", labels(1), types(1), act, con), ...
            side("B", labels(2), types(2), act, con)];
row.AtoB_nConnections = between(con, labels(1), labels(2));
row.BtoA_nConnections = between(con, labels(2), labels(1));
row.syncModalPairCount = S.observedGlobalMaximumPairCount;
row.syncNull95 = S.nullGlobalMaximum95;
row.syncP = S.globalEmpiricalPValue;
row.syncSignificant = double(S.isSignificantAt95);
row.analysedOn = string(datetime('now', 'Format', 'uuuu-MM-dd HH:mm'));
end


function t = side(prefix, label, type, act, con)
% One organoid's half of the row.
k = find(act.organoidLabel == label, 1);
assert(~isempty(k), 'connectoid:NoOrganoid', ...
    'The activity table has no row for %s.', label);
t = table(label, type, act.nUnits(k), act.meanUnitFiringRateHz(k), ...
    act.meanPeakFiringRateHz(k), act.peakFrequencyHz(k), ...
    act.peakToMeanRateRatio(k), between(con, label, label), ...
    'VariableNames', prefix + ["_organoid", "_type", "_nUnits", ...
    "_meanUnitFiringRateHz", "_meanPeakFiringRateHz", "_peakFrequencyHz", ...
    "_peakToMeanRateRatio", "_nWithinConnections"]);
end


function n = between(con, source, target)
k = con.sourceOrganoid == source & con.targetOrganoid == target;
n = sum(con.nConnections(k));
end
