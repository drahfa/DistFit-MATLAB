function m = mixture_biexponential()
%MIXTURE_BIEXPONENTIAL Two-component exponential mixture: w, mu1, mu2
m.name = 'Bi-Exponential';
m.type = 'cont';
m.nParams = 3;
m.paramNames = {'w','mu1','mu2'};
m.support = [0, inf];
m.start = @(x) [0.5, mean(x,'omitnan')/2, mean(x,'omitnan')];
m.fit = @fitf;
m.pdf = @(x,th) th(1).*exppdf_local(x, th(2)) + (1-th(1)).*exppdf_local(x, th(3));
m.cdf = @(x,th) th(1).*expcdf_local(x, th(2)) + (1-th(1)).*expcdf_local(x, th(3));
m.rnd = @(n,th) [rand(n,1) < th(1)].*exprnd(th(2),n,1) + [rand(n,1) >= th(1)].*exprnd(th(3),n,1);

    function th = fitf(x, th0)
        x = x(:); x = x(x>=0);
        if nargin<2 || isempty(th0), th0 = m.start(x); end
        % Bounds: 0<w<1, mu1>0, mu2>0
        nll = @(p) -sum(log(max(m.pdf(x, [min(max(p(1),1e-6),1-1e-6), max(p(2),eps), max(p(3),eps)]), realmin)));
        opts = optimset('display','off');
        p = fminsearch(nll, th0, opts);
        p(1) = min(max(p(1),1e-6),1-1e-6); p(2) = max(p(2),eps); p(3) = max(p(3),eps);
        th = p(:).';
    end
end

function y = exppdf_local(x,mu)
y = zeros(size(x)); mask = x>=0; y(mask) = (1./mu).*exp(-x(mask)./mu);
end
function y = expcdf_local(x,mu)
y = zeros(size(x)); mask = x>=0; y(mask) = 1 - exp(-x(mask)./mu);
end

