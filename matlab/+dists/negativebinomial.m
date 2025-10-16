function m = negativebinomial()
%NEGATIVEBINOMIAL Negative Binomial (r>0, p in (0,1))
m.name = 'Negative Binomial';
m.type = 'disc';
m.nParams = 2;
m.paramNames = {'r','p'};
m.support = [0, inf];
m.start = @(x) start_nbin(x);
m.fit = @fitf;
m.pdf = @(k,th) nbinpdf_local(k, th(1), th(2));
m.cdf = @(k,th) nbincdf_local(k, th(1), th(2));
m.rnd = @(n,th) nbinrnd(max(th(1),eps), min(max(th(2),eps),1-eps), n, 1);

    function th = fitf(x, th0)
        x = x(:); x = x(isfinite(x) & x>=0);
        if isempty(x)
            th = [1, 0.5];
            return;
        end
        if nargin<2 || isempty(th0), th0 = m.start(x); end
        th0 = enforce_bounds(th0);
        phi0 = theta_to_phi(th0);
        opts = optimset('display','off','tolx',1e-8,'tolfun',1e-8,'maxiter',4000,'maxfunevals',8000);
        [phi, ~, exitflag] = fminsearch(@(ph) nll_phi(ph, x), phi0, opts);
        if exitflag <= 0
            % Retry from a variance-inflated guess if needed.
            alt = enforce_bounds([max(th0(1)*2, eps), min(max((th0(2)+0.5)/2, eps), 1-eps)]);
            phi = theta_to_phi(alt);
            [phi, ~] = fminsearch(@(ph) nll_phi(ph, x), phi, opts);
        end
        th = phi_to_theta(phi);
    end

    function phi = theta_to_phi(theta)
        r = max(theta(1), eps);
        p = min(max(theta(2), eps), 1-eps);
        phi = [log(r); log(p/(1-p))];
    end

    function theta = phi_to_theta(phi)
        r = max(exp(phi(1)), eps);
        p = 1./(1+exp(-phi(2)));
        p = min(max(p, eps), 1-eps);
        theta = [r, p];
    end

    function val = nll_phi(phi, data)
        th = phi_to_theta(phi);
        val = -sum(log(max(nbinpdf_local(data, th(1), th(2)), realmin)));
    end

    function out = enforce_bounds(theta)
        out = [max(theta(1), eps), min(max(theta(2), eps), 1-eps)];
    end
end

function th0 = start_nbin(x)
m = mean(x,'omitnan'); v = var(x,'omitnan');
if v <= m
    r = 10; p = min(max(m/(m+v+eps), eps), 1-eps);
else
    p = m/v; p = min(max(p, eps), 1-eps); r = m*p/(1-p);
end
th0 = [r, p];
end

function p = nbinpdf_local(k,r,pp)
kk = floor(k);
p = zeros(size(kk));
mask = kk>=0;
kk = kk(mask);
p(mask) = exp(gammaln(kk+r) - gammaln(kk+1) - gammaln(r) + r*log(pp) + kk.*log(1-pp));
end

function F = nbincdf_local(k,r,pp)
kk = floor(k);
F = betainc(pp, r, kk+1);
end
