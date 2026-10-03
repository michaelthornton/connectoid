function save_figure_set(figures, outputDirectory, names)
% Write each figure as PNG, vector EPS and .fig.
%
% Helvetica and the painters renderer are set before printing, so the EPS
% carries live text Illustrator can edit.
arguments
    figures
    outputDirectory (1,1) string
    names (:,1) string
end
assert(numel(figures)==numel(names),'connectoid:FigureNames', ...
    'Got %d figures and %d names.',numel(figures),numel(names));
for k = 1:numel(figures)
    fig = figures(k);
    assert(isgraphics(fig,'figure'),'connectoid:DeletedFigure', ...
        'Figure "%s" was deleted before it could be exported.',names(k));
    set(fig,'Renderer','painters');
    set(findall(fig,'-property','FontName'),'FontName','Helvetica');
    drawnow;
    stem = char(fullfile(outputDirectory,names(k)));
    print(fig, stem, '-dpng','-r300');
    print(fig, stem, '-depsc2','-vector','-loose');
    savefig(fig, [stem '.fig']);
end
end
