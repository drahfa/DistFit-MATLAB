function results = fitAllDistributions(data, opts)
%FITALLDISTRIBUTIONS Fit multiple distributions and compute GOF metrics.
% results: struct array with fields name, theta, KS, AD, CS, loglik, AIC, BIC
if nargin < 2, opts = struct; end
if ~isfield(opts,'orderBy'), opts.orderBy = 'KS'; end
if ~isfield(opts,'computeP'), opts.computeP = 'none'; end % 'none'|'asymptotic'|'bootstrap'
if ~isfield(opts,'nBootstrap'), opts.nBootstrap = 0; end

x = data(:);
x = x(isfinite(x));
if isempty(x)
    results = struct('name',{},'theta',{},'KS',{},'AD',{},'CS',{},'loglik',{},'AIC',{},'BIC',{});
    return;
end

models = dists.registry();
% If the data are integers, restrict to discrete models
if all(abs(x - round(x)) < 1e-12)
    models = models(strcmp({models.type}, 'disc'));
else
    models = models(strcmp({models.type}, 'cont') | strcmp({models.type}, 'disc'));
end
res = struct('name',{},'theta',{},'KS',{},'AD',{},'CS',{},'KSp',{},'ADp',{},'loglik',{},'AIC',{},'BIC',{});
for k = 1:numel(models)
    M = models(k);
    % Domain check
    if strcmp(M.type,'cont') && M.support(1) >= 0 && any(x < 0)
        % e.g., Exponential/Weibull/Lognormal require x>=0
        xfit = x(x>=0);
        if isempty(xfit), continue; end
    else
        xfit = x;
    end
    % Fit
    try
        th0 = M.start(xfit);
        theta = M.fit(xfit, th0);
    catch
        warning('Fit failed for %s', M.name);
        continue;
    end
    % GOF
    try
        KS = ks_stat(xfit, M.cdf, theta);
    catch, KS = NaN; end
    try
        AD = ad_stat(xfit, M.cdf, theta);
    catch, AD = NaN; end
    % p-values
    KSp = NaN; ADp = NaN;
    switch lower(opts.computeP)
        case 'asymptotic'
            try, KSp = ks_pvalue_asymptotic(KS, numel(xfit)); end
        case 'bootstrap'
            if isfield(M,'rnd') && ~isempty(M.rnd) && opts.nBootstrap>0
                try
                    KSp = bootstrap_pvalue(xfit, @(y) ks_stat(y,M.cdf,theta), M.rnd, theta, opts.nBootstrap);
                    ADp = bootstrap_pvalue(xfit, @(y) ad_stat(y,M.cdf,theta), M.rnd, theta, opts.nBootstrap);
                catch
                end
            end
    end
    try
        CS = chisq_stat(xfit, M.pdf, M.cdf, theta, M.type);
    catch, CS = NaN; end
    % Likelihood and ICs (approx.)
    try
        ll = sum(log(max(M.pdf(xfit, theta), realmin)));
        kpar = numel(theta);
        AIC = 2*kpar - 2*ll;
        BIC = kpar*log(numel(xfit)) - 2*ll;
    catch
        ll = NaN; AIC = NaN; BIC = NaN;
    end
    res(end+1) = struct('name',M.name,'theta',theta,'KS',KS,'AD',AD,'CS',CS, ...
                        'KSp',KSp,'ADp',ADp,'loglik',ll,'AIC',AIC,'BIC',BIC); %#ok<AGROW>
end

% Rank
if isempty(res)
    results = res; return;
end
key = upper(string(opts.orderBy));
vals = arrayfun(@(r) r.(char(key)), res);
[~,idx] = sort(vals, 'ascend', 'MissingPlacement','last');
results = res(idx);
end
