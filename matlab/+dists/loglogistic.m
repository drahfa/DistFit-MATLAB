function m = loglogistic()
%LOGLOGISTIC Log-logistic (Fisk) distribution (alpha>0 scale, beta>0 shape)
m.name = 'Log-Logistic';
m.type = 'cont';
m.nParams = 2;
m.paramNames = {'alpha','beta'};
m.support = [0, inf];
m.start = @(x) [median(x(x>0)) , 1];
m.fit = @fitf;
m.pdf = @(x,th) fisk_pdf(x, th(1), th(2));
m.cdf = @(x,th) fisk_cdf(x, th(1), th(2));
m.rnd = @(n,th) th(1) .* (rand(n,1)./(1-rand(n,1))).^(1/th(2));

    function th = fitf(x, th0)
        x = x(:); x = x(x>0);
        if nargin<2 || isempty(th0), th0 = m.start(x); end
        nll = @(p) -sum(log(max(fisk_pdf(x, max(p(1),eps), max(p(2),eps)), realmin)));
        opts = optimset('display','off');
        p = fminsearch(nll, th0, opts);
        p = max(p, eps);
        th = p(:).';
    end
end

function y = fisk_pdf(x,a,b)
y = zeros(size(x));
mask = x>0 & a>0 & b>0;
xx = x(mask)./a;
y(mask) = (b./a).*xx.^(b-1) ./ (1 + xx.^b).^2;
end

function F = fisk_cdf(x,a,b)
F = zeros(size(x));
mask = x>0 & a>0 & b>0;
xx = x(mask)./a;
F(mask) = 1./(1 + xx.^(-b));
end

