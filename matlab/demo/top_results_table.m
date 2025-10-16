function T = top_results_table(results, k)
%TOP_RESULTS_TABLE Convert top-k fit results to a MATLAB table
if nargin < 2 || isempty(k), k = min(5, numel(results)); end
k = min(k, numel(results));
if k == 0
    T = table; return;
end
names = {results(1:k).name}.';
KS = [results(1:k).KS].';
AD = [results(1:k).AD].';
ChiSq = [results(1:k).CS].';
loglik = [results(1:k).loglik].';
AIC = [results(1:k).AIC].';
BIC = [results(1:k).BIC].';
T = table(names, KS, AD, ChiSq, loglik, AIC, BIC, ...
    'VariableNames', {'Name','KS','AD','ChiSq','LogLik','AIC','BIC'});
end

