function m = lognormal()
%LOGNORMAL Lognormal distribution with mu, sigma of log(X)
m.name = 'Lognormal';
m.type = 'cont';
m.nParams = 2;
m.paramNames = {'mu','sigma'}; % parameters of underlying normal
m.support = [0, inf];
m.start = @(x) start_logn(x);
m.fit = @fitf;
m.pdf = @(x,th) lognpdf_local(x, th(1), th(2));
m.cdf = @(x,th) logncdf_local(x, th(1), th(2));
m.rnd = @(n,th) exp(th(1) + th(2).*randn(n,1));

    function th = fitf(x, th0)
        x = x(:);
        x = x(x>0);
        if nargin < 2 || isempty(th0), th0 = m.start(x); end
        try
            pd = fitdist(x,'Lognormal');
            th = [pd.mu, pd.sigma];
        catch
            th = mle(x, 'distribution','logn', 'start', th0, 'lowerbound',[ -Inf eps ]);
        end
        th = th(:).';
    end
end

function th0 = start_logn(x)
x = x(:);
x = x(x>0);
lx = log(x);
th0 = [mean(lx,'omitnan'), std(lx,0,'omitnan')+eps];
end

function y = lognpdf_local(x,mu,sigma)
y = zeros(size(x));
mask = x>0;
z = (log(x(mask)) - mu)./sigma;
y(mask) = (1./(x(mask).*sigma*sqrt(2*pi))) .* exp(-0.5*z.^2);
end

function y = logncdf_local(x,mu,sigma)
y = zeros(size(x));
mask = x>0;
z = (log(x(mask)) - mu)./(sigma*sqrt(2));
y(mask) = 0.5*(1+erf(z));
end
