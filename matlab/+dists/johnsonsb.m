function m = johnsonsb()
%JOHNSONSB Johnson SB distribution (bounded in (xi, xi+lambda))
m.name = 'Johnson SB';
m.type = 'cont';
m.nParams = 4; % gamma, delta>0, xi, lambda>0
m.paramNames = {'gamma','delta','xi','lambda'};
m.support = [-inf, inf];
m.start = @(x) start_sb(x);
m.fit = @fitf;
m.pdf = @(x,th) sb_pdf(x, th(1), th(2), th(3), th(4));
m.cdf = @(x,th) sb_cdf(x, th(1), th(2), th(3), th(4));
m.rnd = @(n,th) sb_rnd(n, th);

    function th = fitf(x, th0)
        x = x(:); x = x(isfinite(x));
        if nargin<2 || isempty(th0), th0 = m.start(x); end
        % Optimize unconstrained params: [gamma, log(delta), xi, log(lambda)]
        totheta = @(p) [p(1), exp(p(2)), p(3), exp(p(4))];
        nll = @(p) local_nll(p, x, totheta);
        p0 = [th0(1), log(max(th0(2),eps)), th0(3), log(max(th0(4),eps))];
        opts = optimset('display','off');
        p = fminsearch(nll, p0, opts);
        th = totheta(p);
    end
end

function th0 = start_sb(x)
x = x(:); x = x(isfinite(x));
a = min(x); b = max(x); if a==b, b=a+1; end
xi = a - 0.05*(b-a); lambda = (b - a)*1.1; % widen slightly
gamma = 0; delta = 1;
th0 = [gamma, delta, xi, lambda];
end

function f = sb_pdf(x,gamma,delta,xi,lambda)
f = zeros(size(x));
mask = (x>xi) & (x<xi+lambda) & delta>0 & lambda>0;
y = (x(mask) - xi)./lambda;
z = gamma + delta .* log(y./(1-y));
phi = (1/sqrt(2*pi)) * exp(-0.5*z.^2);
f(mask) = phi .* (delta ./ (lambda .* y .* (1-y)));
end

function F = sb_cdf(x,gamma,delta,xi,lambda)
F = zeros(size(x));
mask = (x>xi) & (x<xi+lambda) & delta>0 & lambda>0;
y = (x(mask) - xi)./lambda;
z = gamma + delta .* log(y./(1-y));
F(mask) = 0.5*(1+erf(z./sqrt(2)));
F(x<=xi) = 0; F(x>=xi+lambda) = 1;
end

function r = sb_rnd(n, th)
gamma=th(1); delta=th(2); xi=th(3); lambda=th(4);
Z = randn(n,1);
W = (Z - gamma)./delta;
Y = 1./(1 + exp(-W));
r = xi + lambda.*Y;
end
function val = local_nll(p, x, totheta)
th = totheta(p);
val = -sum(log(max(sb_pdf(x, th(1), th(2), th(3), th(4)), realmin)));
end

