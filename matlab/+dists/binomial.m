function m = binomial()
%BINOMIAL Binomial distribution (n, p)
m.name = 'Binomial';
m.type = 'disc';
m.nParams = 2;
m.paramNames = {'n','p'};
m.support = [0, inf];
m.start = @(x) start_binom(x);
m.fit = @fitf;
m.pdf = @(k,th) binopdf_local(k, th(1), th(2));
m.cdf = @(k,th) binocdf_local(k, th(1), th(2));
m.rnd = @(n,th) binornd(max(round(th(1)),1), min(max(th(2),eps),1-eps), n, 1);

    function th = fitf(x, th0)
        x = x(:); x = x(isfinite(x) & x>=0);
        if nargin<2 || isempty(th0), th0 = m.start(x); end
        % grid search over n
        nmax = max(max(x), ceil(mean(x)*5)+1);
        nmax = max(nmax, 1);
        best = Inf; bestn = th0(1); bestp = th0(2);
        for n = max(1, max(x)):nmax
            p = mean(x)/n; p = min(max(p, eps), 1-eps);
            nll = -sum(log(max(binopdf_local(x, n, p), realmin)));
            if nll < best, best = nll; bestn = n; bestp = p; end
        end
        th = [bestn, bestp];
    end
end

function th0 = start_binom(x)
m = mean(x,'omitnan'); v = var(x,'omitnan');
if v < m*(1-1e-3)
    p = max(eps, 1 - v/max(m,eps)); n = max(1, round(m/p));
else
    n = max(1, round(max(x))); p = min(max(m/n, eps), 1-eps);
end
th0 = [n, p];
end

function p = binopdf_local(k,n,pp)
kk = floor(k);
p = zeros(size(kk));
mask = kk>=0 & kk<=n;
if any(mask)
    kk = kk(mask);
    p(mask) = exp(gammaln(n+1) - gammaln(kk+1) - gammaln(n-kk+1) + kk.*log(pp) + (n-kk).*log(1-pp));
end
end

function F = binocdf_local(k,n,pp)
kk = floor(k);
F = zeros(size(kk));
for i=1:numel(kk)
    if kk(i) < 0
        F(i) = 0;
    elseif kk(i) >= n
        F(i) = 1;
    else
        vals = 0:kk(i);
        F(i) = sum(exp(gammaln(n+1) - gammaln(vals+1) - gammaln(n-vals+1) + vals.*log(pp) + (n-vals).*log(1-pp)));
    end
end
end

