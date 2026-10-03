function name = recording_label(file)
% The name a recording goes into a table under: 000123/data.raw.h5 is
% "000123". A file not called data.raw.h5 is named by its own stem.
arguments
    file (1,1) string
end
[folder, stem] = fileparts(char(file));
if startsWith(stem, 'data.raw')
    [~, name] = fileparts(folder);
else
    name = stem;
end
name = string(name);
end
