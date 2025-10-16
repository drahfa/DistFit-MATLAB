function m = chisquare()
%CHISQUARE Chi-square distribution (df k>0)
m.name = 'Chi-Squared';
m.type = 'cont';
m.nParams = 1;
m.paramNames = {'k'};
m.support = [0, inf];
m.start = @(x) max(eps, mean(x,'omitnan'));
m.fit = @fitf;
m.pdf = @(x,th) gampdf_local(x, th(1)/2, 2);
m.cdf = @(x,th) gamcdf_local(x, th(1)/2, 2);
m.rnd = @(n,th) chi2rnd(max(th(1),eps), n, 1);

    function th = fitf(x, th0)
        x = x(:); x = x(x>=0);
        if nargin<2 || isempty(th0), th0 = m.start(x); end
        % Optimize k only
        nll = @(k) -sum(log(max(gampdf_local(x, max(k,eps)/2, 2), realmin)));
        opts = optimset('display','off');
        k = fminsearch(nll, th0, opts);
        th = max(k, eps);
    end
end

function y = gampdf_local(x,k,theta)
y = zeros(size(x));
mask = x>=0 & k>0 & theta>0;
y(mask) = (x(mask).^(k-1) .* exp(-x(mask)./theta)) ./ (gamma(k).*theta.^k);
end

function F = gamcdf_local(x,k,theta)
F = zeros(size(x));
mask = x>=0 & k>0 & theta>0;
F(mask) = gammainc(x(mask)./theta, k,'lower');
end

