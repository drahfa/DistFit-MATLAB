function m = geometric()
%GEOMETRIC Geometric distribution (failures before first success), p
m.name = 'Geometric';
m.type = 'disc';
m.nParams = 1;
m.paramNames = {'p'};
m.support = [0, inf];
m.start = @(x) 1/(mean(x,'omitnan')+1);
m.fit = @fitf;
m.pdf = @(k,th) geopdf_local(k, th(1));
m.cdf = @(k,th) geocdf_local(k, th(1));
m.rnd = @(n,th) geornd(min(max(th(1),eps),1-eps), n, 1);

    function th = fitf(x, th0)
        x = x(:); x = x(x>=0);
        if nargin<2 || isempty(th0), th0 = m.start(x); end
        p = 1/(mean(x,'omitnan')+1); p = min(max(p, eps), 1-eps);
        th = p(:).';
    end
end

function p = geopdf_local(k,pp)
kk = floor(k); p = zeros(size(kk)); mask = kk>=0; p(mask) = pp.*(1-pp).^kk(mask);
end
function F = geocdf_local(k,pp)
kk = floor(k); F = zeros(size(kk)); mask = kk>=0; F(mask) = 1 - (1-pp).^(kk(mask)+1);
end

