function m = gumbel()
%GUMBEL Extreme Value (Type I, maxima) aka Gumbel
m.name = 'Gumbel';
m.type = 'cont';
m.nParams = 2;
m.paramNames = {'mu','sigma'};
m.support = [-inf, inf];
m.start = @(x) [mean(x,'omitnan') - 0.57721*std(x,0,'omitnan'), std(x,0,'omitnan')*sqrt(6)/pi];
m.fit = @fitf;
m.pdf = @(x,th) evpdf_local(x, th(1), th(2));
m.cdf = @(x,th) evcdf_local(x, th(1), th(2));
m.rnd = @(n,th) th(1) - th(2).*log(-log(rand(n,1)));

    function th = fitf(x, th0)
        x = x(:);
        if nargin<2 || isempty(th0), th0 = m.start(x); end
        try
            [mu, sigma] = evfit(x);
            th = [mu, sigma];
        catch
            nll = @(p) -sum(log(max(evpdf_local(x, p(1), max(p(2),eps)), realmin)));
            opts = optimset('display','off');
            p = fminsearch(nll, th0, opts);
            p(2) = max(p(2), eps);
            th = p(:).';
        end
    end
end

function y = evpdf_local(x,mu,sigma)
z = (x-mu)./sigma;
e = exp(-z);
y = (1./sigma).*e.*exp(-e);
end

function F = evcdf_local(x,mu,sigma)
z = (x-mu)./sigma;
F = exp(-exp(-z));
end

