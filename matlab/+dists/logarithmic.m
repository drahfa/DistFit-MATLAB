function m = logarithmic()
%LOGARITHMIC Log-series distribution (p in (0,1))
m.name = 'Logarithmic';
m.type = 'disc';
m.nParams = 1;
m.paramNames = {'p'};
m.support = [1, inf];
m.start = @(x) min(max(1 - exp(-mean(x,'omitnan')), eps), 1-eps);
m.fit = @fitf;
m.pdf = @(k,th) logser_pdf(k, th(1));
m.cdf = @(k,th) logser_cdf(k, th(1));
m.rnd = @(n,th) logser_rnd(n, th(1));

    function th = fitf(x, th0)
        x = x(:); x = x(x>=1);
        if nargin<2 || isempty(th0), th0 = m.start(x); end
        nll = @(p) -sum(log(max(logser_pdf(x, min(max(p,eps),1-eps)), realmin)));
        opts = optimset('display','off');
        p = fminsearch(nll, th0, opts);
        th = min(max(p, eps), 1-eps);
    end
end

function p = logser_pdf(k,pp)
kk = floor(k);
p = zeros(size(kk));
mask = kk>=1 & pp>0 & pp<1;
kkm = kk(mask);
Z = -1/log(1-pp);
p(mask) = Z * (pp.^kkm) ./ kkm;
end

function F = logser_cdf(k,pp)
kk = floor(k); F = zeros(size(kk));
for i=1:numel(kk)
    if kk(i) < 1
        F(i) = 0;
    else
        vals = 1:kk(i);
        F(i) = sum(logser_pdf(vals, pp));
        F(i) = min(F(i),1);
    end
end
end

function r = logser_rnd(n,pp)
r = zeros(n,1);
u = rand(n,1);
% Inverse CDF via cumulative sum until exceed u
for i=1:n
    k=1; c=0; 
    while true
        c = c + logser_pdf(k, pp);
        if u(i) <= c || k>1e5, r(i)=k; break; end
        k = k + 1;
    end
end
end

