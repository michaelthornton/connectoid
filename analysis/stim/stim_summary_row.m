function row = stim_summary_row(label, recs, R, M, opts)
% One connectoid as one row: both directions and the call for each.
%
% Column names carry no organoid names, so one table holds any pair; stimA and
% stimB say which was which. stimA stimulated the first organoid in the ROI
% table. Every value is read from the organoid that was not stimulated.
arguments
    label (1,1) string
    recs (1,2) struct
    R table
    M (1,2) struct
    opts.WindowSec (1,2) double = [0.030 0.100]
end
types = R.type(:)';
driven = strings(1, 2);
for p = 1:2
    [~, d] = stim_organoids(recs(p), R);
    driven(p) = d.type;
end
assert(numel(unique(driven)) == 2, 'stim:SameDriven', ...
    'Both recordings stimulated %s; they are not the two directions.', driven(1));
a = find(driven == types(1), 1);
b = 3 - a;

row = table(label, types(1), types(2), ...
    'VariableNames', {'dataset', 'stimA', 'stimB'});
row = [row, side("stimA", recs(a), M(a)), side("stimB", recs(b), M(b))];
row.outcome = outcomeOf(M(a).ciSeparated, M(b).ciSeparated, types);
row.windowMs = sprintf('%g-%g', 1000 * opts.WindowSec(1), 1000 * opts.WindowSec(2));
row.analysedOn = string(datetime('now', 'Format', 'uuuu-MM-dd HH:mm'));
end


function t = side(prefix, S, m)
% One stimulation direction. connected is 0 or 1 rather than a logical, so a
% reread table concatenates with a fresh row.
t = table(S.recording, string(S.file), numel(S.onsets), m.respondingOrganoid, ...
    m.nElectrodes, m.baselineSpS, m.deltaPct, m.deltaPctCI, m.pDelta, ...
    double(m.ciSeparated), m.fracAbove, m.pFrac, ...
    'VariableNames', prefix + ["_recording", "_file", "_nTrials", ...
        "_respondingOrganoid", "_nElectrodes", "_baselineHz", ...
        "_windowPctOfBaseline", "_ci95", "_pDelta", "_connected", ...
        "_pctTrialsResponding", "_pFrac"]);
end


function s = outcomeOf(a, b, types)
% The four exclusive outcome categories.
if a && b
    s = types(1) + "<->" + types(2);
elseif a
    s = types(1) + "->" + types(2) + " only";
elseif b
    s = types(2) + "->" + types(1) + " only";
else
    s = "no response";
end
end
