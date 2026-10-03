function T = append_summary_row(file, row)
% Add one recording to a running table. A row whose dataset name is already in
% the file is replaced, so reanalysing a recording leaves one row.
arguments
    file (1,1) string
    row table
end
if isfile(file)
    % Dataset names like 000123 are text; readtable would make them 123.
    isText = varfun(@(v) isstring(v) || ischar(v) || iscellstr(v), row, ...
        'OutputFormat', 'uniform');
    io = detectImportOptions(file, 'TextType', 'string');
    text = intersect(io.VariableNames, row.Properties.VariableNames(isText));
    if ~isempty(text), io = setvartype(io, text, 'string'); end
    T = readtable(file, io);

    old = string(T.Properties.VariableNames);
    new = string(row.Properties.VariableNames);
    assert(isequal(old, new), 'connectoid:ColumnMismatch', ...
        ['%s has different columns from this run, so the rows do not belong ' ...
         'in one table. Only in it: %s. Only in this run: %s.'], ...
        file, join(setdiff(old, new), ", "), join(setdiff(new, old), ", "));
    T(T.dataset == row.dataset, :) = [];
    T = [T; row];
else
    T = row;
end
T = sortrows(T, 'dataset');
writetable(T, file);
fprintf('  %s now holds %d recording(s)\n', file, height(T));
end
