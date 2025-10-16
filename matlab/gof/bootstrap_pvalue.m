function p = bootstrap_pvalue(x, statfun, rndfun, theta, nboot)
%BOOTSTRAP_PVALUE Parametric bootstrap p-value for a statistic
% x: data, statfun: @(y)->stat, rndfun: @(n,theta)->y, theta: fitted params
% nboot: number of resamples
if nargin<5 || isempty(nboot) || nboot<=0, p = NaN; return; end
x = x(:); x = x(isfinite(x));
if isempty(x) || isempty(rndfun), p = NaN; return; end
obs = statfun(x);
cnt = 0; n = numel(x);
for b=1:nboot
    y = rndfun(n, theta);
    s = statfun(y);
    if s >= obs, cnt = cnt + 1; end
end
p = (cnt + 1) / (nboot + 1);
end

