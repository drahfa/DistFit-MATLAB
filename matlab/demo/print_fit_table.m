function print_fit_table(titleStr, results, k)
%PRINT_FIT_TABLE Display a formatted table of fit metrics
if nargin < 3 || isempty(k), k = min(5, numel(results)); end
k = min(k, numel(results));
fprintf('\n=== %s ===\n', titleStr);
fprintf('%-4s  %-22s  %8s  %8s  %8s\n','Rank','Name','KS','AD','ChiSq');
for i = 1:k
    r = results(i);
    fprintf('#%-3d  %-22s  %8.4g  %8.4g  %8.4g\n', i, r.name, r.KS, r.AD, r.CS);
end
end

