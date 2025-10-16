function qqplot_model(x, M, theta)
%QQPLOT_MODEL Q-Q plot using numeric inverse CDF if needed
x = x(:); x = x(isfinite(x));
x = sort(x);
n = numel(x);
p = ((1:n)' - 0.5)/n;
% Inverse CDF
if isfield(M,'inv') && ~isempty(M.inv)
    q = M.inv(p, theta);
else
    q = inv_cdf_numeric(p, @(t) M.cdf(t,theta), min(x), max(x));
end
plot(q, x, '.'); hold on;
mn = min([x; q]); mx = max([x; q]);
plot([mn mx],[mn mx],'r--');
xlabel('Theoretical Quantiles'); ylabel('Data Quantiles');
title(sprintf('%s Q-Q plot', M.name));
hold off;
end

function q = inv_cdf_numeric(p, cdfh, xmin, xmax)
% Numerically invert CDF on [xmin-5*range, xmax+5*range]
if ~isfinite(xmin) || ~isfinite(xmax) || xmin==xmax
    xmin = 0; xmax = 1; 
end
range = xmax - xmin;
a = xmin - 5*range; b = xmax + 5*range;
q = zeros(size(p));
for i=1:numel(p)
    f = @(t) cdfh(t) - p(i);
    try
        q(i) = fzero(f, [a b]);
    catch
        % fallback coarse search
        xs = linspace(a,b,1000); [~,idx] = min(abs(cdfh(xs)-p(i))); q(i) = xs(idx);
    end
end
end

