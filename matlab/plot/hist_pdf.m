function hist_pdf(x, M, theta, nbins)
%HIST_PDF Plot histogram and overlay model PDF
if nargin<4 || isempty(nbins), nbins = max(20, ceil(sqrt(numel(x)))); end
x = x(:); x = x(isfinite(x));
histogram(x, nbins, 'Normalization','pdf'); hold on;
xmin = min(x); xmax = max(x);
xx = linspace(xmin, xmax, 200);
yy = M.pdf(xx, theta);
plot(xx, yy, 'r-', 'LineWidth', 2);
title(sprintf('%s PDF fit', M.name));
hold off;
end

