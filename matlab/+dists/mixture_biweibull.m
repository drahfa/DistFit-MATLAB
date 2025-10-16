function m = mixture_biweibull()
%MIXTURE_BIWEIBULL Two-component Weibull mixture: w, a1, b1, a2, b2
m.name = 'Bi-Weibull';
m.type = 'cont';
m.nParams = 5;
m.paramNames = {'w','A1','B1','A2','B2'};
m.support = [0, inf];
m.start = @(x) start_params(x);
m.fit = @fitf;
m.pdf = @(x,th) th(1).*wblpdf_local(x, th(2), th(3)) + (1-th(1)).*wblpdf_local(x, th(4), th(5));
m.cdf = @(x,th) th(1).*wblcdf_local(x, th(2), th(3)) + (1-th(1)).*wblcdf_local(x, th(4), th(5));
m.rnd = @(n,th) ([rand(n,1)<th(1)].*wblrnd_local(th(2),th(3),n,1) + [rand(n,1)>=th(1)].*wblrnd_local(th(4),th(5),n,1));

    function th = fitf(x, th0)
        x = x(:); x = x(x>=0);
        if nargin<2 || isempty(th0), th0 = m.start(x); end
        nll = @(p) -sum(log(max(m.pdf(x, [min(max(p(1),1e-6),1-1e-6), max(p(2),eps), max(p(3),eps), max(p(4),eps), max(p(5),eps)]), realmin)));
        opts = optimset('display','off');
        p = fminsearch(nll, th0, opts);
        p(1) = min(max(p(1),1e-6),1-1e-6);
        p(2:5) = max(p(2:5), eps);
        th = p(:).';
    end
end

function th0 = start_params(x)
% Seed from single Weibull fit
try
    pd = fitdist(x,'Weibull'); a=pd.A; b=pd.B;
catch
    a = mean(x); b = 1.5;
end
w = 0.5;
th0 = [w, a*0.7, max(b*0.8,0.5), a*1.3, max(b*1.2,0.5)];
end

function y = wblpdf_local(x,a,b)
y = zeros(size(x)); mask = x>=0; y(mask) = (b./a).*(x(mask)./a).^(b-1).*exp(-(x(mask)./a).^b);
end
function y = wblcdf_local(x,a,b)
y = zeros(size(x)); mask = x>=0; y(mask) = 1 - exp(-(x(mask)./a).^b);
end
function r = wblrnd_local(a,b,m,n)
U = rand(m,n); r = a .* (-log(max(1e-12,1-U))).^(1./b);
end

