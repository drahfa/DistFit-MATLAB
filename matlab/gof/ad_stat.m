function A2 = ad_stat(x, cdfh, theta)
%AD_STAT Anderson-Darling statistic (no small-sample correction)
x = x(:);
x = x(isfinite(x));
if isempty(x); A2 = NaN; return; end
x = sort(x);
n = numel(x);
F = cdfh(x, theta);
% Clamp to avoid log(0)
F = min(max(F, eps), 1-eps);
i = (1:n)';
A2 = -n - mean((2*i-1).*(log(F) + log(1 - F(end:-1:1))));
end

