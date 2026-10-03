function S = stim_read_recording(file, opts)
% Spike table, routed map, stimulation sites and train onsets from a MaxWell
% stimulation recording. Takes the H5 or the folder holding it.
%
% Reads no voltage traces, so the MaxWell HDF5 filter plugin is not needed.
%
%   S = stim_read_recording("/path/to/000124")
arguments
    file (1,1) string
    opts.WellID (1,1) double {mustBeInteger,mustBePositive} = 1
end
file = resolve_recording(file);   % the H5, or the folder holding it
G = sprintf('/data_store/data%04d/', opts.WellID - 1);
f = char(file);
try
    h5info(f, G);
catch
    error('stim:NoWell', 'Well %d is not in this file (no %s).\n%s', ...
        opts.WellID, G, f);
end

m = h5read(f, [G 'settings/mapping']);
fsAll = h5read(f, [G 'settings/sampling']);
fs = double(fsAll(1));
f0 = double(h5read(f, [G 'groups/routed/frame_nos'], 1, 1));
frameInfo = h5info(f, [G 'groups/routed/frame_nos']);
nFrames = double(frameInfo.Dataspace.Size);
sp = h5read(f, [G 'spikes']);
ev = h5read(f, [G 'events']);
raw = h5read(f, '/assay/inputs/electrodes');
assay = jsondecode(raw{1});

% Stimulation sites appear more than once in settings/mapping, so dedupe or
% their spikes are counted twice.
rows = unique([double(m.channel(:)) double(m.electrode(:))], 'rows');
ch2el = -ones(max(rows(:,1)) + 2, 1);
ch2el(rows(:,1) + 1) = rows(:,2);

% One row per unique electrode, with its routed position.
[uel, ia] = unique(double(m.electrode(:)));
S.file = file;
S.recording = recording_label(file);
S.electrodes = uel;
S.x = double(m.x(ia));
S.y = double(m.y(ia));

% MaxOne is a 220-by-120 grid on a 17.5 um pitch, so the electrode id fixes
% the position. Checked rather than assumed.
assert(max(abs(mod(uel, 220) * 17.5 - S.x)) < 1e-6 && ...
       max(abs(floor(uel / 220) * 17.5 - S.y)) < 1e-6, ...
    'stim:Geometry', '%s: the id-to-position formula disagrees with the routed map.', ...
    S.recording);

S.fs = fs;
S.duration = nFrames / fs;
S.t = (double(sp.frameno) - f0) / fs;
S.el = ch2el(double(sp.channel) + 1);
S.amp = abs(double(sp.amplitude));

fn = fieldnames(assay.electrodes);
assert(numel(fn) >= opts.WellID, 'stim:NoAssayWell', ...
    '%s lists %d well(s) in the assay; well %d was asked for.', ...
    S.recording, numel(fn), opts.WellID);
well = assay.electrodes.(fn{opts.WellID});
assert(isstruct(well) && isfield(well, 'stimulation'), 'stim:NotStimFile', ...
    '%s has no stimulation electrode list; it is not a stimulation recording.', ...
    S.recording);
S.stimElectrodes = double(well.stimulation(:));

types = double(ev.eventtype);
pulses = sort((double(ev.frameno(types == 1)) - f0) / fs);
assert(~isempty(pulses), 'stim:NoEvents', '%s has no stimulation events.', S.recording);
% Grouped by time, not eventid: some acquisitions reuse ids across trials.
starts = [1; find(diff(pulses) > 0.015) + 1];
S.onsets = pulses(starts);
S.lastPulse = pulses([starts(2:end) - 1; numel(pulses)]);
S.pulses = pulses;
end
