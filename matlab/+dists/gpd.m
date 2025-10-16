function m = gpd()
%GPD Generalized Pareto distribution (sigma>0, xi)
m.name = 'Generalized Pareto';
m.type = 'cont';
m.nParams = 2;
m.paramNames = {'sigma','xi'};
m.support = [0, inf];
m.start = @(x) start_gpd(x);
m.fit = @fitf;
m.pdf = @(x,th) gpd_pdf(x, th(1), th(2));
m.cdf = @(x,th) gpd_cdf(x, th(1), th(2));
m.rnd = @(n,th) gpd_rnd(n, th(1), th(2));

    function th = fitf(x, th0)
        x = x(:); x = x(isfinite(x));
        x = x(x >= 0);
        if isempty(x)
            th = [1, 0];
            return;
        end
        if nargin < 2 || isempty(th0), th0 = m.start(x); end
        th0 = enforce_bounds(th0);
        try %#ok<TRYNC>
            [khat, sigmahat] = gpfit(x);
            th = enforce_bounds([sigmahat, khat]);
            return;
        end
        phi0 = theta_to_phi(th0);
        opts = optimset('display','off','tolx',1e-8,'tolfun',1e-8,'maxiter',4000,'maxfunevals',8000);
        [phi, fval, exitflag] = fminsearch(@(ph) nll_phi(ph, x), phi0, opts);
        if exitflag <= 0 || ~isfinite(fval)
            alt = enforce_bounds([th0(1)*1.2, th0(2)/2]);
            [phi, ~] = fminsearch(@(ph) nll_phi(ph, x), theta_to_phi(alt), opts);
        end
        th = enforce_bounds(phi_to_theta(phi));
    end

    function phi = theta_to_phi(theta)
        sigma = max(theta(1), eps);
        xi = clamp(theta(2), -0.5 + 1e-6, 0.5 - 1e-6);
        phi = [log(sigma); atanh(xi / 0.5)];
    end

    function theta = phi_to_theta(phi)
        sigma = max(exp(phi(1)), eps);
        xi = clamp(0.5 * tanh(phi(2)), -0.5 + 1e-6, 0.5 - 1e-6);
        theta = [sigma, xi];
    end

    function val = nll_phi(phi, data)
        th = phi_to_theta(phi);
        logpdf = gpd_logpdf(data, th(1), th(2));
        if any(~isfinite(logpdf))
            val = Inf;
        else
            val = -sum(logpdf);
        end
    end

    function th = enforce_bounds(theta)
        sigma = max(theta(1), eps);
        xi = clamp(theta(2), -0.5 + 1e-6, 0.5 - 1e-6);
        th = [sigma, xi];
    end
end

function th0 = start_gpd(x)
x = x(:); x = x(isfinite(x));
x = x(x >= 0);
if isempty(x)
    th0 = [1, 0];
    return;
end
sigma = std(x);
if ~isfinite(sigma) || sigma <= 0
    sigma = median(abs(x - median(x)))*1.4826;
    if ~isfinite(sigma) || sigma <= 0, sigma = 1; end
end
sigma = max(sigma, eps);
xi = skewness(x);
if ~isfinite(xi), xi = 0; end
xi = clamp(0.1 * xi, -0.3, 0.3);
th0 = [sigma, xi];
end

function logpdf = gpd_logpdf(x, sigma, xi)
sigma = max(sigma, eps);
x = x(:);
if abs(xi) < 1e-8
    logpdf = -log(sigma) - x ./ sigma;
    logpdf(x < 0) = -Inf;
    return;
end
z = 1 + xi .* (x ./ sigma);
logpdf = -Inf(size(x));
mask = (x >= 0) & (z > 0);
logpdf(mask) = -log(sigma) - (1./xi + 1).*log(z(mask));
end

function y = gpd_pdf(x, sigma, xi)
sigma = max(sigma, eps);
y = zeros(size(x));
mask = x >= 0;
if ~any(mask)
    return;
end
if abs(xi) < 1e-8
    y(mask) = (1./sigma) .* exp(-x(mask) ./ sigma);
    return;
end
z = 1 + xi .* (x(mask) ./ sigma);
valid = z > 0;
indices = find(mask);
y(indices(valid)) = (1./sigma) .* z(valid).^(-1./xi - 1);
end

function F = gpd_cdf(x, sigma, xi)
sigma = max(sigma, eps);
F = zeros(size(x));
mask = x >= 0;
if ~any(mask)
    return;
end
if abs(xi) < 1e-8
    F(mask) = 1 - exp(-x(mask) ./ sigma);
    return;
end
z = 1 + xi .* (x(mask) ./ sigma);
valid = z > 0;
indices = find(mask);
F(indices(valid)) = 1 - z(valid).^(-1./xi);
F(indices(~valid)) = 1;
end

function r = gpd_rnd(n, sigma, xi)
sigma = max(sigma, eps);
if numel(n) == 1
    U = rand(n,1);
else
    U = rand(n);
end
U = max(min(U, 1 - realmin), realmin);
if abs(xi) < 1e-8
    r = -sigma .* log(1 - U);
else
    r = sigma/xi * ((1 - U).^(-xi) - 1);
end
end

function val = clamp(x, lo, hi)
val = min(max(x, lo), hi);
end

