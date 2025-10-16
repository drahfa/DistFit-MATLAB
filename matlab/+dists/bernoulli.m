function m = bernoulli()
%BERNOULLI Bernoulli distribution (p)
m.name = 'Bernoulli';
m.type = 'disc';
m.nParams = 1;
m.paramNames = {'p'};
m.support = [0, 1];
m.start = @(x) min(max(mean(x,'omitnan'),eps),1-eps);
m.fit = @fitf;
m.pdf = @(k,th) (k==0).*(1-th(1)) + (k==1).*th(1);
m.cdf = @(k,th) (k<0).*0 + (k<1).*(1-th(1)) + (k>=1).*1;
m.rnd = @(n,th) rand(n,1) < th(1);

    function th = fitf(x, th0)
        x = x(:); x = x(x==0 | x==1);
        if nargin<2 || isempty(th0), th0 = m.start(x); end
        th = min(max(mean(x,'omitnan'),eps),1-eps);
        th = th(:).';
    end
end

