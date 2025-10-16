function m = studentt()
%STUDENTT t location-scale (mu, sigma>0, nu>0)
m.name = "Student's t";
m.type = 'cont';
m.nParams = 3;
m.paramNames = {'mu','sigma','nu'};
m.support = [-inf, inf];
m.start = @(x) [mean(x,'omitnan'), std(x,0,'omitnan')+eps, 5];
m.fit = @fitf;
m.pdf = @(x,th) tpdf((x-th(1))./th(2), th(3))./th(2);
m.cdf = @(x,th) tcdf((x-th(1))./th(2), th(3));
m.rnd = @(n,th) th(1) + th(2).*trnd(th(3), n, 1);

    function th = fitf(x, th0)
        x = x(:);
        if nargin<2 || isempty(th0), th0 = m.start(x); end
        try
            pd = fitdist(x, 'tLocationScale');
            th = [pd.mu, pd.sigma, pd.nu];
        catch
            ll = @(p) -sum(log(max(tpdf((x-p(1))./max(p(2),eps), max(p(3),eps))./max(p(2),eps), realmin)));
            opts = optimset('display','off');
            p = fminsearch(ll, th0, opts);
            p(2) = max(p(2), eps); p(3) = max(p(3), eps);
            th = p(:).';
        end
    end
end

