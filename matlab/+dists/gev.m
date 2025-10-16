function m = gev()
%GEV Generalized Extreme Value distribution (mu, sigma>0, kappa)
m.name = 'GEV';
m.type = 'cont';
m.nParams = 3;
m.paramNames = {'mu','sigma','kappa'};
m.support = [-inf, inf];
m.start = @(x) start_gev(x);
m.fit = @fitf;
m.pdf = @(x,th) gevpdf_local(x, th(1), th(2), th(3));
m.cdf = @(x,th) gevcdf_local(x, th(1), th(2), th(3));
m.rnd = @(n,th) gevrnd_local(n, th(1), th(2), th(3));

    function th = fitf(x, th0)
        x = x(:); x = x(isfinite(x));
        if isempty(x)
            th = [0, 1, 0];
            return;
        end
        if nargin < 2 || isempty(th0), th0 = m.start(x); end
        th0 = enforce_bounds(th0);
        try %#ok<TRYNC>
            [khat, sigmahat, muhat] = gevfit(x);
            th = enforce_bounds([muhat, sigmahat, khat]);
            return;
        end
        phi0 = theta_to_phi(th0);
        opts = optimset('display','off','tolx',1e-8,'tolfun',1e-8,'maxiter',4000,'maxfunevals',8000);
        [phi, fval, exitflag] = fminsearch(@(ph) nll_phi(ph, x), phi0, opts);
        if exitflag <= 0 || ~isfinite(fval)
            alt = enforce_bounds([th0(1), th0(2)*1.3, th0(3)/2]);
            [phi, ~] = fminsearch(@(ph) nll_phi(ph, x), theta_to_phi(alt), opts);
        end
        th = enforce_bounds(phi_to_theta(phi));
    end

    function phi = theta_to_phi(theta)
        mu = theta(1);
        sigma = max(theta(2), eps);
        kappa = clamp(theta(3), -0.5 + 1e-6, 0.5 - 1e-6);
        ratio = max(min(kappa / 0.5, 0.999), -0.999);
        phi = [mu; log(sigma); atanh(ratio)];
    end

    function theta = phi_to_theta(phi)
        mu = phi(1);
        sigma = max(exp(phi(2)), eps);
        kappa = clamp(0.5 * tanh(phi(3)), -0.5 + 1e-6, 0.5 - 1e-6);
        theta = [mu, sigma, kappa];
    end

    function val = nll_phi(phi, data)
        th = phi_to_theta(phi);
        logpdf = gev_logpdf(data, th(1), th(2), th(3));
        if any(~isfinite(logpdf))
            val = Inf;
        else
            val = -sum(logpdf);
        end
    end

    function th = enforce_bounds(theta)
        mu = theta(1);
        sigma = max(theta(2), eps);
        kappa = clamp(theta(3), -0.5 + 1e-6, 0.5 - 1e-6);
        th = [mu, sigma, kappa];
    end
end

function th0 = start_gev(x)
x = x(:); x = x(isfinite(x));
if isempty(x)
    th0 = [0, 1, 0];
    return;
end
mu = mean(x);
if ~isfinite(mu)
    mu = median(x);
end
sigma = std(x);
if ~isfinite(sigma) || sigma <= 0
    x0 = x - mu;
    mad0 = median(abs(x0));
    sigma = max(mad0 * 1.4826, 1);
end
sigma = max(sigma, eps);
z = (x - mu) ./ sigma;
sk = mean(z.^3);
if ~isfinite(sk)
    sk = 0;
end
kappa = clamp(0.15*sign(sk), -0.3, 0.3);
th0 = [mu, sigma, kappa];
end

function logpdf = gev_logpdf(x, mu, sigma, kappa)
sigma = max(sigma, eps);
z = (x - mu) ./ sigma;
t = 1 + kappa .* z;
if abs(kappa) < 1e-6
    logpdf = -log(sigma) - z - exp(-z);
    return;
end
if any(t <= 0)
    logpdf = -Inf(size(x));
    return;
end
logt = log(t);
logpdf = -log(sigma) - (1./kappa + 1).*logt - t.^(-1./kappa);
end

function y = gevpdf_local(x, mu, sigma, kappa)
sigma = max(sigma, eps);
z = (x - mu) ./ sigma;
t = 1 + kappa .* z;
y = zeros(size(x));
mask = t > 0;
if abs(kappa) < 1e-6
    y(mask) = (1./sigma) .* exp(-z(mask) - exp(-z(mask)));
else
    y(mask) = (1./sigma) .* t(mask).^(-1./kappa - 1) .* exp(-t(mask).^(-1./kappa));
end
end

function F = gevcdf_local(x, mu, sigma, kappa)
sigma = max(sigma, eps);
z = (x - mu) ./ sigma;
t = 1 + kappa .* z;
F = zeros(size(x));
if abs(kappa) < 1e-6
    F = exp(-exp(-z));
    return;
end
mask = t > 0;
F(mask) = exp(-t(mask).^(-1./kappa));
if kappa > 0
    F(~mask) = 0;
else
    F(~mask) = 1;
end
end

function r = gevrnd_local(n, mu, sigma, kappa)
sigma = max(sigma, eps);
if numel(n) == 1
    U = rand(n, 1);
else
    U = rand(n);
end
U = max(min(U, 1 - realmin), realmin);
if abs(kappa) < 1e-6
    r = mu - sigma .* log(-log(U));
else
    r = mu + (sigma./kappa) .* ((-log(U)).^(-kappa) - 1);
end
end

function val = clamp(x, lo, hi)
val = min(max(x, lo), hi);
end
