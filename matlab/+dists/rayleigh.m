function m = rayleigh()
%RAYLEIGH Rayleigh distribution (sigma>0)
m.name = 'Rayleigh';
m.type = 'cont';
m.nParams = 1;
m.paramNames = {'sigma'};
m.support = [0, inf];
m.start = @(x) sqrt(mean(x(x>=0).^2,'omitnan')/2);
m.fit = @fitf;
m.pdf = @(x,th) raylpdf_local(x, th(1));
m.cdf = @(x,th) raylcdf_local(x, th(1));
m.rnd = @(n,th) th(1).*sqrt(-2*log(max(1e-12,1-rand(n,1))));

    function th = fitf(x, th0)
        x = x(:); x = x(x>=0);
        if nargin<2 || isempty(th0), th0 = m.start(x); end
        try
            pd = fitdist(x,'Rayleigh');
            th = pd.B; % sigma
        catch
            th = mle(x, 'pdf', @(y,s) raylpdf_local(y,s), 'start', th0, 'lowerbound', eps);
        end
        th = th(:).';
    end
end

function y = raylpdf_local(x,s)
y = zeros(size(x));
mask = x>=0 & s>0;
y(mask) = (x(mask)./s.^2).*exp(-(x(mask).^2)./(2*s.^2));
end

function F = raylcdf_local(x,s)
F = zeros(size(x));
mask = x>=0 & s>0;
F(mask) = 1 - exp(-(x(mask).^2)./(2*s.^2));
end

