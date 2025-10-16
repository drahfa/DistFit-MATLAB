function m = weibull()
%WEIBULL Weibull distribution with scale A>0 and shape B>0
m.name = 'Weibull';
m.type = 'cont';
m.nParams = 2;
m.paramNames = {'A','B'}; % MATLAB uses A(scale), B(shape)
m.support = [0, inf];
m.start = @(x) start_wbl(x);
m.fit = @fitf;
m.pdf = @(x,th) wblpdf_local(x, th(1), th(2));
m.cdf = @(x,th) wblcdf_local(x, th(1), th(2));
m.rnd = @(n,th) th(1) .* (-log(max(1e-12,1-rand(n,1)))).^(1./th(2));

    function th = fitf(x, th0)
        x = x(:);
        x = x(x>=0);
        if nargin < 2 || isempty(th0), th0 = m.start(x); end
        th0 = max(th0, eps);
        try
            pd = fitdist(x,'Weibull');
            th = [pd.A, pd.B];
        catch
            th = mle(x, 'distribution','wbl', 'start', th0, 'lowerbound',[eps eps]);
        end
        th = th(:).';
    end
end

function th0 = start_wbl(x)
x = x(:);
x = x(x>0);
if isempty(x)
    th0 = [1,1];
    return;
end
% Method-of-moments rough init
m1 = mean(x,'omitnan');
m2 = var(x,'omitnan');
cv2 = m2/(m1^2 + eps);
% Approximate shape from CV using simple mapping
k = max(0.1, min(10, 1/(cv2+eps))); % crude
lambda = m1 / gamma(1+1/k);
th0 = [lambda, k];
end

function y = wblpdf_local(x,a,b)
y = zeros(size(x));
mask = x>=0;
y(mask) = (b./a) .* (x(mask)./a).^(b-1) .* exp(-(x(mask)./a).^b);
end

function y = wblcdf_local(x,a,b)
y = zeros(size(x));
mask = x>=0;
y(mask) = 1 - exp(-(x(mask)./a).^b);
end
