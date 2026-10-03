function t = stim_norm_text(ax, x, y, str, varargin)
% Place text in normalized axes units.
%
% text(ax, x, y, str, 'Units', 'normalized') does not do this: Position is set
% from x and y before Units, so the coordinates are read as data units. Create
% first, then set Units and Position.
t = text(ax, 0, 0, str, varargin{:});
set(t, 'Units', 'normalized', 'Position', [x y 0]);
end
