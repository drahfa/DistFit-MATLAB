function m = inversegaussian()
%INVERSEGAUSSIAN Inverse Gaussian (Wald) distribution (mu>0, lambda>0)
m.name = 'Inverse Gaussian';
m.type = 'cont';
m.nParams = 2;
m.paramNames = {'mu','lambda'};
m.support = [0, inf];
m.start = @(x) [mean(x(x>0),'omitnan'), 1/var(1./x(x>0),'omitnan')];
m.fit = @fitf;
m.pdf = @(x,th) ig_pdf(x, th(1), th(2));
m.cdf = @(x,th) ig_cdf(x, th(1), th(2));
m.rnd = [];

    function th = fitf(x, th0)
        x = x(:); x = x(x>0);
        if nargin<2 || isempty(th0), th0 = m.start(x); end
        try
            pd = fitdist(x, 'InverseGaussian');
            th = [pd.mu, pd.lambda];
        catch
            nll = @(p) -sum(log(max(ig_pdf(x, max(p(1),eps), max(p(2),eps)), realmin)));
            opts = optimset('display','off');
            p = fminsearch(nll, th0, opts);
            p = max(p, eps);
            th = p(:).';
        end
    end
end

function y = ig_pdf(x,mu,lambda)
y = zeros(size(x));
mask = x>0 & mu>0 & lambda>0;
xx = x(mask);
y(mask) = sqrt(lambda ./ (2*pi*xx.^3)) .* exp(-(lambda.*(xx-mu).^2)./(2*mu.^2.*xx));
end

function F = ig_cdf(x,mu,lambda)
F = zeros(size(x));
mask = x>0 & mu>0 & lambda>0;
xx = x(mask);
z1 = sqrt(lambda./xx) .* ((xx./mu) - 1);
z2 = -sqrt(lambda./xx) .* ((xx./mu) + 1);
F(mask) = 0.5*(erfc(-z1./sqrt(2)) + exp(2*lambda./mu) .* erfc(-z2./sqrt(2)));
F(mask) = min(max(F(mask),0),1);
end
