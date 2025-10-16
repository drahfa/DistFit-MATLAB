function m = hypergeometric()
%HYPERGEOMETRIC Hypergeometric distribution (N,K,n)
m.name = 'Hypergeometric';
m.type = 'disc';
m.nParams = 3;
m.paramNames = {'N','K','n'};
m.support = [0, inf];
m.start = @(x) start_hyper(x);
m.fit = @fitf;
m.pdf = @(k,th) hygepdf_local(k, round(th(1)), round(th(2)), round(th(3)));
m.cdf = @(k,th) hygecdf_local(k, round(th(1)), round(th(2)), round(th(3)));
m.rnd = @(n,th) hygernd(round(th(1)), round(th(2)), round(th(3)), n, 1);

    function th = fitf(x, th0)
        x = x(:); x = x(x>=0);
        if nargin<2 || isempty(th0), th0 = m.start(x); end
        % Grid search around starts (coarse), clamp validity
        N0 = max(round(th0(1)), max(x)+1);
        n0 = max(round(th0(3)), max(x));
        K0 = max(round(th0(2)), 1);
        best = Inf; bestTheta = [N0,K0,n0];
        Ngrid = unique(max(max(x)+1, N0 + (-5:5)));
        ngrid = unique(max(max(x), n0 + (-5:5)));
        Kgrid = unique(max(1, K0 + (-5:5)));
        for N=Ngrid
            for n=ngrid
                if n>N, continue; end
                for K=Kgrid
                    if K>N, continue; end
                    if any(x > min(K,n)), continue; end
                    p = hygepdf_local(x, N, K, n); p(p<=0)=realmin;
                    nll = -sum(log(p));
                    if nll < best, best=nll; bestTheta=[N,K,n]; end
                end
            end
        end
        th = bestTheta;
    end
end

function th0 = start_hyper(x)
m = mean(x,'omitnan'); v = var(x,'omitnan');
% Crude heuristics: n ~ max(x), N ~ 10*n, K from mean relation
n = max(1, round(max(x)));
N = max(n+1, 10*n);
K = min(N-1, max(1, round(m * N / max(n,1))));
th0 = [N, K, n];
end

function p = hygepdf_local(x,N,K,n)
k = floor(x);
p = zeros(size(k));
mask = k>=0 & k<=min(K,n) & N>=K & N>=n;
kk = k(mask);
p(mask) = exp(gammaln(K+1)-gammaln(kk+1)-gammaln(K-kk+1) + ...
               gammaln(N-K+1)-gammaln(n-kk+1)-gammaln(N-K-(n-kk)+1) - ...
               (gammaln(N+1)-gammaln(n+1)-gammaln(N-n+1)));
end

function F = hygecdf_local(x,N,K,n)
kk = floor(x);
F = zeros(size(kk));
for i=1:numel(kk)
    if kk(i) < 0
        F(i) = 0;
    else
        vals = 0:min(kk(i),min(K,n));
        F(i) = sum(hygepdf_local(vals, N, K, n));
        F(i) = min(F(i),1);
    end
end
end

