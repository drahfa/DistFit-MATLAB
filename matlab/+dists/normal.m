function m = normal()
%NORMAL Normal distribution model definition
m.name = 'Normal';
m.type = 'cont';
m.nParams = 2;
m.paramNames = {'mu','sigma'};
m.support = [-inf, inf];
m.start = @(x) [mean(x,'omitnan'), std(x,0,'omitnan')+eps];
m.fit = @fitf;
m.pdf = @(x,th) normpdf_local(x, th(1), th(2));
m.cdf = @(x,th) normcdf_local(x, th(1), th(2));
m.rnd = @(n,th) th(1) + th(2).*randn(n,1);

    function th = fitf(x, th0)
        x = x(:);
        if nargin < 2 || isempty(th0), th0 = m.start(x); end
        % Try fitdist, fallback to mle
        try
            pd = fitdist(x, 'Normal');
            th = [pd.mu, pd.sigma];
        catch
            th = mle(x, 'distribution','normal', 'start', th0);
        end
    end
end

function y = normpdf_local(x,mu,sigma)
z = (x-mu)./sigma;
y = (1./(sqrt(2*pi).*sigma)) .* exp(-0.5*z.^2);
end

function y = normcdf_local(x,mu,sigma)
z = (x-mu)./(sigma*sqrt(2));
y = 0.5*(1+erf(z));
end
