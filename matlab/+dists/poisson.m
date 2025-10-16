function m = poisson()
%POISSON Poisson distribution for nonnegative integers
m.name = 'Poisson';
m.type = 'disc';
m.nParams = 1;
m.paramNames = {'lambda'};
m.support = [0, inf];
m.start = @(x) max(eps, mean(x,'omitnan'));
m.fit = @fitf;
m.pdf = @(x,th) poisspdf_local(x, th(1));
m.cdf = @(x,th) poisscdf_local(x, th(1));
m.rnd = @(n,th) poissrnd(max(th(1),eps), n, 1);

    function th = fitf(x, th0)
        x = x(:);
        x = x(x>=0 & isfinite(x));
        if nargin < 2 || isempty(th0), th0 = m.start(x); end
        try
            pd = fitdist(x, 'Poisson');
            th = pd.lambda;
        catch
            % MLE for Poisson: lambda = mean(x)
            th = mean(x,'omitnan');
        end
        th = max(th, eps);
        th = th(:).';
    end
end

function p = poisspdf_local(k,lambda)
k = floor(k);
p = zeros(size(k));
mask = k>=0 & isfinite(k);
kk = k(mask);
p(mask) = exp(kk.*log(lambda) - lambda - gammaln(kk+1));
end

function F = poisscdf_local(k,lambda)
k = floor(k);
F = zeros(size(k));
mask = k>=0 & isfinite(k);
kk = k(mask);
% CDF via incomplete gamma relation or summation
maxk = max(kk(:));
vals = 0:maxk;
pmf = exp(vals.*log(lambda) - lambda - gammaln(vals+1));
cum = cumsum(pmf);
F(mask) = cum(kk+1);
end
