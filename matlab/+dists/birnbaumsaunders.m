function m = birnbaumsaunders()
%BIRNBAUMSAUNDERS Fatigue Life (alpha>0, beta>0)
m.name = 'Birnbaum-Saunders';
m.type = 'cont';
m.nParams = 2;
m.paramNames = {'alpha','beta'};
m.support = [0, inf];
m.start = @(x) [0.5, mean(x(x>0),'omitnan')];
m.fit = @fitf;
m.pdf = @(x,th) bs_pdf(x, th(1), th(2));
m.cdf = @(x,th) bs_cdf(x, th(1), th(2));
m.rnd = [];

    function th = fitf(x, th0)
        x = x(:); x = x(x>0);
        if nargin<2 || isempty(th0), th0 = m.start(x); end
        try
            pd = fitdist(x, 'BirnbaumSaunders');
            th = [pd.alpha, pd.beta];
        catch
            nll = @(p) -sum(log(max(bs_pdf(x, max(p(1),eps), max(p(2),eps)), realmin)));
            opts = optimset('display','off');
            p = fminsearch(nll, th0, opts);
            p = max(p, eps);
            th = p(:).';
        end
    end
end

function y = bs_pdf(x,alpha,beta)
y = zeros(size(x));
mask = x>0 & alpha>0 & beta>0;
xx = x(mask);
z = (sqrt(xx./beta) - sqrt(beta./xx))./alpha;
y(mask) = (1./(2*alpha*sqrt(2*pi)*xx)) .* (sqrt(xx./beta) + sqrt(beta./xx)) .* exp(-0.5*z.^2);
end

function F = bs_cdf(x,alpha,beta)
F = zeros(size(x));
mask = x>0 & alpha>0 & beta>0;
xx = x(mask);
z = (sqrt(xx./beta) - sqrt(beta./xx))./alpha;
F(mask) = 0.5*(1+erf(z./sqrt(2)));
end
