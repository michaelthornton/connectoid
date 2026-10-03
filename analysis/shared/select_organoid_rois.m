function positions = select_organoid_rois(electrodeData, labels, useSquareRois, ...
    descriptions)
% Draw, edit, and confirm one organoid box per label.
% Returns one [x y width height] row per label.
%
% descriptions is what to call each organoid on screen, for callers whose labels
% are not meaningful on their own. It defaults to the labels.

arguments
    electrodeData table
    labels (:, 1) string
    useSquareRois (1, 1) logical = true
    descriptions (:, 1) string = labels
end
if numel(descriptions) ~= numel(labels)
    error('select_organoid_rois:DescriptionCount', ...
        'Got %d descriptions for %d organoid labels.', ...
        numel(descriptions), numel(labels));
end

requiredVariables = {'x', 'y', 'nSpikes'};
missingVariables = setdiff(requiredVariables, electrodeData.Properties.VariableNames);
if ~isempty(missingVariables)
    error('select_organoid_rois:InvalidElectrodeTable', ...
        'Electrode table is missing: %s.', strjoin(missingVariables, ', '));
end
if isempty(labels)
    error('select_organoid_rois:MissingLabels', ...
        'At least one organoid label is required.');
end

fig = figure('Color', 'w', 'Name', 'Define organoid ROIs', ...
    'NumberTitle', 'off', 'WindowStyle', 'normal', ...
    'Position', [180 120 1000 700]);
ax = axes(fig, 'Position', [0.08 0.12 0.88 0.76]);
hold(ax, 'on');
inactive = electrodeData.nSpikes == 0;
scatter(ax, electrodeData.x(inactive), electrodeData.y(inactive), 8, ...
    [0.85 0.85 0.85], 'filled');
scatter(ax, electrodeData.x(~inactive), electrodeData.y(~inactive), 14, ...
    'k', 'filled');
axis(ax, 'equal');
axis(ax, 'ij');
axis(ax, 'tight');
xlabel(ax, 'x (\mum)');
ylabel(ax, 'y (\mum)');

colors = lines(max(numel(labels), 1));
instruction = uicontrol(fig, 'Style', 'text', 'Units', 'normalized', ...
    'Position', [0.08 0.91 0.88 0.07], 'BackgroundColor', 'w', ...
    'HorizontalAlignment', 'left', 'FontSize', 12, 'FontWeight', 'bold');
confirmButton = uicontrol(fig, 'Style', 'pushbutton', 'String', 'Confirm ROIs', ...
    'Units', 'normalized', 'Position', [0.70 0.025 0.13 0.055], ...
    'FontWeight', 'bold', 'Enable', 'off', ...
    'Callback', @(~, ~) confirmSelection(fig));
cancelButton = uicontrol(fig, 'Style', 'pushbutton', 'String', 'Cancel', ...
    'Units', 'normalized', 'Position', [0.84 0.025 0.10 0.055], ...
    'Enable', 'off', 'Callback', @(~, ~) cancelSelection(fig));
set(fig, 'CloseRequestFcn', @(~, ~) cancelSelection(fig));
setappdata(fig, 'RoiSelectionConfirmed', false);
% Bring it to the front: it is a blocking dialog, and behind the command
% window it looks like MATLAB has hung.
movegui(fig, 'center');
figure(fig);
drawnow;
assignment = strjoin(compose("Organoid %d = %s", (1:numel(labels))', descriptions), '; ');
fprintf('ROI assignment: %s.\n', assignment);

roiHandles = gobjects(numel(labels), 1);
for iOrganoid = 1:numel(labels)
    instruction.String = sprintf( ...
        ['%s. Step %d of %d: drag a box around %s. ', ...
        'Release the mouse to continue.'], ...
        assignment, iOrganoid, numel(labels), descriptions(iOrganoid));
    drawnow;
    roiHandles(iOrganoid) = drawrectangle(ax, 'FixedAspectRatio', useSquareRois, ...
        'Label', descriptions(iOrganoid), 'LineWidth', 2, ...
        'Color', colors(iOrganoid, :));
end

instruction.String = sprintf(['%s. Every box is editable. ', ...
    'Drag or resize them, then click Confirm ROIs.'], assignment);
confirmButton.Enable = 'on';
cancelButton.Enable = 'on';
title(ax, 'Organoid ROI assignment');
drawnow;
uiwait(fig);

if ~isvalid(fig) || ~getappdata(fig, 'RoiSelectionConfirmed')
    if isvalid(fig)
        delete(fig);
    end
    error('select_organoid_rois:SelectionCancelled', ...
        'Organoid ROI selection was cancelled.');
end

positions = zeros(numel(labels), 4);
for iOrganoid = 1:numel(labels)
    if ~isvalid(roiHandles(iOrganoid))
        delete(fig);
        error('select_organoid_rois:MissingRoi', ...
            'The %s ROI was deleted before confirmation.', descriptions(iOrganoid));
    end
    positions(iOrganoid, :) = double(roiHandles(iOrganoid).Position(:).');
end
delete(fig);
end

function confirmSelection(fig)
if isvalid(fig)
    setappdata(fig, 'RoiSelectionConfirmed', true);
    uiresume(fig);
end
end

function cancelSelection(fig)
if isvalid(fig)
    setappdata(fig, 'RoiSelectionConfirmed', false);
    uiresume(fig);
end
end
