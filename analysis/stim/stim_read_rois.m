function R = stim_read_rois(file)
% The two organoid rectangles and their types, from a CSV written by a
% previous run.
%
% Columns: type, x, y, width, height, one row per organoid, in micrometres.
% The type column is required; row order carries no identity, since the boxes
% are drawn in whatever order they were clicked.
arguments
    file (1,1) string
end
assert(isfile(file), 'stim:NoRoiFile', 'No ROI file at:\n%s', file);
T = readtable(file, 'TextType', 'string');
need = ["type" "x" "y" "width" "height"];
have = string(T.Properties.VariableNames);
assert(all(ismember(need, have)), 'stim:RoiColumns', ...
    '%s needs the columns %s; it has %s.', file, strjoin(need, ", "), strjoin(have, ", "));
assert(height(T) == 2, 'stim:RoiCount', '%s has %d rows; two are expected.', ...
    file, height(T));
assert(~any(ismissing(T.type)) && ~any(strlength(T.type) == 0), 'stim:RoiTypes', ...
    '%s has a blank type; row order does not say which organoid is which.', file);
R = table(string(T.type), T.x, T.y, T.width, T.height, ...
          'VariableNames', {'type','x','y','width','height'});
end
