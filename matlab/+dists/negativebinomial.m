function m = negativebinomial()
%NEGATIVEBINOMIAL Negative Binomial (r>0, p in (0,1))
m.name = 'Negative Binomial';
m.type = 'disc';
m.nParams = 2;
m.paramNames = {'r','p'};
m.support = [0, inf];
m.start = @(x) start_nbin(x);
m.fit = @fitf;
m.pdf = @(k,th) nbinpdf_local(k, th(1), th(2));
m.cdf = @(k,th) nbincdf_local(k, th(1), th(2));
m.rnd = @(n,th) nbinrnd(max(th(1),eps), min(max(th(2),eps),1-eps), n, 1);

    function th = fitf(x, th0)
        x = x(:); x = x(isfinite(x) & x>=0);
        if nargin<2 || isempty(th0), th0 = m.start(x); end
        try
            phat = mle(x, 'pdf', @(y,r,p) nbinpdf_local(y,r,p), 'start', th0, 'lowerbound',[eps eps], 'upperbound',[Inf 1-eps]);
            th = phat(:).';
        catch
            nll = @(p) -sum(log(max(nbinpdf_local(x, max(p(1),eps), min(max(p(2),eps),1-eps)), realmin)));
            opts = optimset('display','off');
            p = fminsearch(nll, th0, opts);
            th = [max(p(1),eps), min(max(p(2),eps),1-eps)];
        end
    end
end

function th0 = start_nbin(x)
m = mean(x,'omitnan'); v = var(x,'omitnan');
if v <= m
    r = 10; p = min(max(m/(m+v+eps), eps), 1-eps);
else
    p = m/v; p = min(max(p, eps), 1-eps); r = m*p/(1-p);
end
th0 = [r, p];
end

function p = nbinpdf_local(k,r,pp)
kk = floor(k);
p = zeros(size(kk));
mask = kk>=0;
kk = kk(mask);
p(mask) = exp(gammaln(kk+r) - gammaln(kk+1) - gammaln(r) + r*log(pp) + kk.*log(1-pp));
end

function F = nbincdf_local(k,r,pp)
kk = floor(k);
F = betainc(pp, r, kk+1);
end

