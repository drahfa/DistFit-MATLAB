function m = johnsonsu()
%JOHNSONSU Johnson SU distribution (unbounded)
m.name = 'Johnson SU';
m.type = 'cont';
m.nParams = 4; % gamma, delta>0, xi, lambda>0
m.paramNames = {'gamma','delta','xi','lambda'};
m.support = [-inf, inf];
m.start = @(x) start_su(x);
m.fit = @fitf;
m.pdf = @(x,th) su_pdf(x, th(1), th(2), th(3), th(4));
m.cdf = @(x,th) su_cdf(x, th(1), th(2), th(3), th(4));
m.rnd = @(n,th) su_rnd(n, th);

    function th = fitf(x, th0)
        x = x(:); x = x(isfinite(x));
        if nargin<2 || isempty(th0), th0 = m.start(x); end
        totheta = @(p) [p(1), exp(p(2)), p(3), exp(p(4))];
        nll = @(p) local_nll_su(p, x, totheta);
        p0 = [th0(1), log(max(th0(2),eps)), th0(3), log(max(th0(4),eps))];
        opts = optimset('display','off');
        p = fminsearch(nll, p0, opts);
        th = totheta(p);
    end
end

function th0 = start_su(x)
mu = median(x,'omitnan'); s = std(x,0,'omitnan');
if ~isfinite(s) || s<=0, s=1; end
gamma = 0; delta = 1; xi = mu; lambda = s;
th0 = [gamma, delta, xi, lambda];
end

function f = su_pdf(x,gamma,delta,xi,lambda)
f = zeros(size(x));
mask = isfinite(x) & delta>0 & lambda>0;
z = gamma + delta .* asinh((x(mask) - xi)./lambda);
phi = (1/sqrt(2*pi)) * exp(-0.5*z.^2);
f(mask) = phi .* (delta ./ sqrt((x(mask)-xi).^2 + lambda.^2));
end

function F = su_cdf(x,gamma,delta,xi,lambda)
F = zeros(size(x));
mask = isfinite(x) & delta>0 & lambda>0;
z = gamma + delta .* asinh((x(mask) - xi)./lambda);
F(mask) = 0.5*(1+erf(z./sqrt(2)));
end

function r = su_rnd(n, th)
gamma=th(1); delta=th(2); xi=th(3); lambda=th(4);
Z = randn(n,1);
W = (Z - gamma)./delta;
r = xi + lambda.*sinh(W);
end
function val = local_nll_su(p, x, totheta)
th = totheta(p);
val = -sum(log(max(su_pdf(x, th(1), th(2), th(3), th(4)), realmin)));
end

