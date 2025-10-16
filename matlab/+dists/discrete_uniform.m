function m = discrete_uniform()
%DISCRETE_UNIFORM Discrete uniform on integers [a,b]
m.name = 'Discrete Uniform';
m.type = 'disc';
m.nParams = 2;
m.paramNames = {'a','b'};
m.support = [-inf, inf];
m.start = @(x) [floor(min(x)), ceil(max(x))];
m.fit = @fitf;
m.pdf = @(k,th) du_pdf(k, round(th(1)), round(th(2)));
m.cdf = @(k,th) du_cdf(k, round(th(1)), round(th(2)));
m.rnd = @(n,th) randi([round(th(1)), round(th(2))], n, 1);

    function th = fitf(x, th0)
        a = floor(min(x)); b = ceil(max(x));
        if a==b, b=a+1; end
        th = [a,b];
    end
end

function p = du_pdf(k,a,b)
kk = floor(k);
width = b-a+1; p = zeros(size(kk));
mask = kk>=a & kk<=b;
p(mask) = 1/width;
end

function F = du_cdf(k,a,b)
kk = floor(k); F = zeros(size(kk));
F(kk<a) = 0; F(kk>=b) = 1;
mask = kk>=a & kk<b;
F(mask) = (kk(mask)-a+1)/(b-a+1);
end

