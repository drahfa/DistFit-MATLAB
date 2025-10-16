function X2 = chisq_stat(x, pdfh, cdfh, theta, type)
%CHISQ_STAT Pearson Chi-square statistic using simple binning.
% For 'cont', use Freedman–Diaconis bins. For 'disc', use unique values.
x = x(:);
x = x(isfinite(x));
if isempty(x); X2 = NaN; return; end

switch type
    case 'disc'
        % Group by integer values
        k = floor(x);
        [cats,~,idx] = unique(k);
        O = accumarray(idx, 1);
        E = numel(x) * pdfh(cats, theta);
        % Avoid zero expected
        E(E<=0) = eps;
        X2 = sum((O - E).^2 ./ E);
    otherwise
        % Continuous: FD rule for bins
        edges = bin_edges_fd(x);
        % Ensure at least 5 bins
        if numel(edges) < 6
            nb = max(5, ceil(sqrt(numel(x))));
            edges = linspace(min(x), max(x), nb+1);
        end
        O = histcounts(x, edges);
        % Expected from CDF difference
        Fe = cdfh(edges, theta);
        Fe = min(max(Fe,0),1);
        Ee = diff(Fe) * numel(x);
        Ee(Ee<=0) = eps;
        X2 = sum((O - Ee).^2 ./ Ee);
end
end

function edges = bin_edges_fd(x)
% Freedman–Diaconis rule
x = x(:);
x = sort(x);
iqr = iqr_local(x);
if iqr <= 0
    edges = unique(x)';
    if numel(edges) < 2
        edges = [x(1)-0.5, x(1)+0.5];
    end
    return;
end
binw = 2*iqr*(numel(x)^(-1/3));
if binw <= 0, binw = (max(x)-min(x))/(ceil(sqrt(numel(x)))+eps); end
edges = min(x):binw:max(x);
if edges(end) < max(x)
    edges = [edges, max(x)];
end
end

function v = iqr_local(x)
q1 = quantile(x, 0.25);
q3 = quantile(x, 0.75);
v = q3 - q1;
end

