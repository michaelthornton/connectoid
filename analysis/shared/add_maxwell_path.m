function toolboxRoot = add_maxwell_path(projectRoot)
%ADD_MAXWELL_PATH Put the MaxWell MATLAB toolbox on the path.
%   Looks for a toolbox already on the path first, so a user who installed it
%   with MaxLab Live does not have to copy it into this project. Falls back to
%   the copy bundled at <projectRoot>/Maxwell/matlab. Returns the folder that
%   was added, or "" when the toolbox was already available.
%
%   The toolbox is MaxWell Biosystems' and is not ours to maintain; see
%   README.md for the version these scripts were written against.

if nargin < 1 || strlength(string(projectRoot)) == 0
    projectRoot = fileparts(fileparts(mfilename('fullpath')));
end

if exist('mxw.fileManager', 'class') == 8
    toolboxRoot = "";
    return
end

bundled = fullfile(projectRoot, 'Maxwell', 'matlab');
if isfolder(bundled)
    addpath(genpath(bundled));
    toolboxRoot = string(bundled);
    return
end

error('connectoid:MissingMaxwellToolbox', ...
    ['The MaxWell MATLAB toolbox was not found.\n' ...
     'Either add your own copy to the MATLAB path (it ships with MaxLab ' ...
     'Live), or place it at:\n  %s\n' ...
     'See README.md for the version these scripts were written against.'], ...
    bundled);
end
