function m = logistic()
%LOGISTIC Logistic distribution (mu, s>0)
m.name = 'Logistic';
m.type = 'cont';
m.nParams = 2;
m.paramNames = {'mu','s'};
m.support = [-inf, inf];
m.start = @(x) start_logistic(x);
m.fit = @fitf;
m.pdf = @(x,th) logistic_pdf(x, th(1), th(2));
m.cdf = @(x,th) logistic_cdf(x, th(1), th(2));
m.rnd = @(n,th) th(1) + th(2).*log(rand(n,1)./(1-rand(n,1)));

    function th = fitf(x, th0)
        x = x(:);
        if nargin<2 || isempty(th0), th0 = m.start(x); end
        ll = @(p) -sum(log(max(logistic_pdf(x,p(1),max(p(2),eps)), realmin)));
        opts = optimset('display','off');
        p = fminsearch(ll, th0, opts);
        p(2) = max(p(2), eps);
        th = p(:).';
    end
end

function th0 = start_logistic(x)
mu = median(x,'omitnan');
sigma = std(x,0,'omitnan');
s = sigma * sqrt(3)/pi;
if ~isfinite(s) || s<=0, s = 1; end
th0 = [mu, s];
end

function y = logistic_pdf(x,mu,s)
z = (x-mu)./s;
ez = exp(-z);
den = s.*(1+ez).^2;
y = ez ./ den;
end

function F = logistic_cdf(x,mu,s)
z = (x-mu)./s;
F = 1./(1+exp(-z));
end

