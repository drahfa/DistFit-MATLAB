function m = exponential()
%EXPONENTIAL Exponential distribution with mean mu > 0
m.name = 'Exponential';
m.type = 'cont';
m.nParams = 1;
m.paramNames = {'mu'};
m.support = [0, inf];
m.start = @(x) mean(x(x>=0),'omitnan');
m.fit = @fitf;
m.pdf = @(x,th) exppdf_local(x, th(1));
m.cdf = @(x,th) expcdf_local(x, th(1));
m.rnd = @(n,th) exprnd(th(1), n, 1);

    function th = fitf(x, th0)
        x = x(:);
        x = x(x>=0);
        if nargin < 2 || isempty(th0) || ~isfinite(th0)
            th0 = m.start(x);
        end
        th0 = max(th0, eps);
        try
            pd = fitdist(x, 'Exponential');
            th = pd.mu;
        catch
            th = mle(x, 'distribution','exp', 'start', th0);
        end
        th = th(:).';
    end
end

function y = exppdf_local(x,mu)
y = zeros(size(x));
mask = x>=0;
y(mask) = (1./mu).*exp(-x(mask)./mu);
end

function y = expcdf_local(x,mu)
y = zeros(size(x));
mask = x>=0;
y(mask) = 1 - exp(-x(mask)./mu);
end
