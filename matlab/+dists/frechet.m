function m = frechet()
%FRECHET Frechet (Type II extreme value) (scale s>0, shape k>0)
m.name = 'Frechet';
m.type = 'cont';
m.nParams = 2;
m.paramNames = {'s','k'};
m.support = [0, inf];
m.start = @(x) [median(x(x>0)), 1];
m.fit = @fitf;
m.pdf = @(x,th) frechet_pdf(x, th(1), th(2));
m.cdf = @(x,th) frechet_cdf(x, th(1), th(2));
m.rnd = @(n,th) th(1) .* (-log(rand(n,1))).^(-1/th(2));

    function th = fitf(x, th0)
        x = x(:); x = x(x>0);
        if nargin<2 || isempty(th0), th0 = m.start(x); end
        nll = @(p) -sum(log(max(frechet_pdf(x, max(p(1),eps), max(p(2),eps)), realmin)));
        opts = optimset('display','off');
        p = fminsearch(nll, th0, opts);
        p = max(p, eps);
        th = p(:).';
    end
end

function y = frechet_pdf(x,s,k)
y = zeros(size(x));
mask = x>0 & s>0 & k>0;
xx = x(mask)./s;
y(mask) = (k./s) .* xx.^(-1-k) .* exp(-xx.^(-k));
end

function F = frechet_cdf(x,s,k)
F = zeros(size(x));
mask = x>0 & s>0 & k>0;
xx = x(mask)./s;
F(mask) = exp(-xx.^(-k));
end

