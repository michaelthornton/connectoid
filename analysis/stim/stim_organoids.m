function [undriven, driven] = stim_organoids(S, R)
% Split the ROIs into the driven and undriven organoid.
%
% Driven = the ROI containing the stimulation sites, which are read from the
% assay. Reported values come from the other organoid: the stimulation sites
% and their neighbours lose most of their detections for the whole file.
%
% Each output is a struct with .type, .electrodes and .nElectrodes.
inRoi = @(k) S.electrodes(S.x >= R.x(k) & S.x <= R.x(k) + R.width(k) & ...
                          S.y >= R.y(k) & S.y <= R.y(k) + R.height(k));
e1 = inRoi(1); e2 = inRoi(2);
has1 = any(ismember(S.stimElectrodes, e1));
has2 = any(ismember(S.stimElectrodes, e2));
assert(xor(has1, has2), 'stim:StimSitesSplit', ...
    ['%s: the stimulation sites fall in %s, which is meant to be exactly one ROI. ' ...
     'The ROIs and the assay disagree; do not report this recording.'], ...
    S.recording, string(has1 * 1 + has2 * 2));
if has1
    driven = mk(R.type(1), e1); undriven = mk(R.type(2), e2);
else
    driven = mk(R.type(2), e2); undriven = mk(R.type(1), e1);
end
end

function o = mk(typ, els)
o = struct('type', typ, 'electrodes', els, 'nElectrodes', numel(els));
assert(o.nElectrodes > 0, 'stim:EmptyRoi', 'No routed electrodes inside the %s ROI.', typ);
end
