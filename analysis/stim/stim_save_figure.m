function stim_save_figure(fig, stem)
% Write a figure in the four formats worth keeping.
%
%   .eps  live text, lossless embedded images. Use this one in Illustrator.
%   .pdf  live text, but exportgraphics writes images as downsampled JPEG.
%   .svg  MATLAB outlines all text, so nothing in it is editable.
%   .png  flat raster for slides.
%
% Resolution also sets the resampling of embedded images in the vector formats.
set(fig, 'Color', 'w', 'InvertHardcopy', 'off');
exportgraphics(fig, stem + ".eps", 'ContentType', 'vector', 'Resolution', 600);
exportgraphics(fig, stem + ".pdf", 'ContentType', 'vector', 'Resolution', 600);
print(fig, char(stem + ".svg"), '-dsvg');
exportgraphics(fig, stem + ".png", 'Resolution', 300);
fprintf('  wrote %s.eps / .pdf / .svg / .png\n', stem);
end
