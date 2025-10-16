function m = gamma()
%GAMMA Gamma distribution (shape k>0, scale theta>0)
m.name = 'Gamma';
m.type = 'cont';
m.nParams = 2;
m.paramNames = {'k','theta'};
m.support = [0, inf];
m.start = @(x) start_gamma(x);
m.fit = @fitf;
m.pdf = @(x,th) gampdf_local(x, th(1), th(2));
m.cdf = @(x,th) gamcdf_local(x, th(1), th(2));
m.rnd = @(n,th) gamrnd(th(1), th(2), n, 1);

    function th = fitf(x, th0)
        x = x(:); x = x(x>=0);
        if nargin<2 || isempty(th0), th0 = m.start(x); end
        try
            pd = fitdist(x,'Gamma'); th = [pd.a, pd.b];
        catch
            th = mle(x, 'pdf', @(y,k,theta) gampdf_local(y,k,theta), ...
                'start', th0, 'lowerbound',[eps eps]);
        end
        th = th(:).';
    end
end

function th0 = start_gamma(x)
x = x(:); x = x(x>0);
mu = mean(x,'omitnan'); v = var(x,'omitnan');
k = max(eps, mu^2/(v+eps)); theta = max(eps, v/(mu+eps));
th0 = [k, theta];
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

