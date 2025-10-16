function D = ks_stat(x, cdfh, theta)
%KS_STAT One-sample KS statistic (no p-value), general CDF handle.
% x: data vector; cdfh: @(x,theta) -> F(x); theta: params
x = x(:);
x = x(isfinite(x));
if isempty(x); D = NaN; return; end
x = sort(x);
n = numel(x);
F = cdfh(x, theta);
F = min(max(F, 0), 1);
% Empirical CDF just after each sample
i = (1:n)';
Femp = i./n;
FempPrev = (i-1)./n;
D = max([max(abs(F - Femp)), max(abs(F - FempPrev))]);
end

