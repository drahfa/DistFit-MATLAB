function m = studentt()
%STUDENTT t location-scale (mu, sigma>0, nu>0)
m.name = "Student's t";
m.type = 'cont';
m.nParams = 3;
m.paramNames = {'mu','sigma','nu'};
m.support = [-inf, inf];
m.start = @(x) [mean(x,'omitnan'), std(x,0,'omitnan')+eps, 5];
m.fit = @fitf;
m.pdf = @(x,th) tpdf((x-th(1))./th(2), th(3))./th(2);
m.cdf = @(x,th) tcdf((x-th(1))./th(2), th(3));
m.rnd = @(n,th) th(1) + th(2).*trnd(th(3), n, 1);

    function th = fitf(x, th0)
        x = x(:); x = x(isfinite(x));
        if isempty(x)
            th = [0, 1, 5];
            return;
        end
        if nargin<2 || isempty(th0), th0 = m.start(x); end
        th0 = enforce_bounds(th0);
        phi0 = theta_to_phi(th0);
        opts = optimset('display','off','tolx',1e-8,'tolfun',1e-8,'maxiter',4000,'maxfunevals',8000);
        [phi, ~, exitflag] = fminsearch(@(ph) nll_phi(ph, x), phi0, opts);
        if exitflag <= 0
            alt = enforce_bounds([th0(1), th0(2)*1.5, max(th0(3)/2, 2)]);
            [phi, ~] = fminsearch(@(ph) nll_phi(ph, x), theta_to_phi(alt), opts);
        end
        th = phi_to_theta(phi);
    end

    function phi = theta_to_phi(theta)
        mu = theta(1);
        sigma = max(theta(2), eps);
        nu = max(theta(3), 1+eps);
        phi = [mu; log(sigma); log(nu-1)];
    end

    function theta = phi_to_theta(phi)
        mu = phi(1);
        sigma = max(exp(phi(2)), eps);
        nu = max(exp(phi(3)) + 1, 1+eps);
        theta = [mu, sigma, nu];
    end

    function val = nll_phi(phi, data)
        th = phi_to_theta(phi);
        z = (data - th(1)) ./ th(2);
        logpdf = log(max(tpdf(z, th(3)) ./ th(2), realmin));
        val = -sum(logpdf);
    end

    function out = enforce_bounds(theta)
        out = [theta(1), max(theta(2), eps), max(theta(3), 1+eps)];
    end
end
