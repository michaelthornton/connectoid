function file = resolve_recording(given)
% The H5 of a recording, given either the file or the folder holding it.
%
% MaxWell writes data.raw.h5 inside a folder named for the recording, so both
% are accepted.
arguments
    given (1,1) string
end
if isfile(given)
    file = given;
    return
end
assert(isfolder(given), 'connectoid:NoRecording', ...
    'No such file or folder:\n%s', given);

inside = fullfile(given, "data.raw.h5");
if isfile(inside)
    file = inside;
    return
end
found = dir(char(fullfile(given, "*.h5")));
found = found(~[found.isdir]);
assert(~isempty(found), 'connectoid:NoRecordingInFolder', ...
    'No .h5 file in:\n%s', given);
assert(isscalar(found), 'connectoid:AmbiguousRecording', ...
    ['%s holds %d .h5 files and no data.raw.h5, so which one to read is not ' ...
     'obvious. Give the file itself.'], given, numel(found));
file = string(fullfile(found.folder, found.name));
end
