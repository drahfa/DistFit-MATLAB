function p = ks_pvalue_asymptotic(D, n)
%KS_PVALUE_ASYMPTOTIC Asymptotic p-value for one-sample KS
if ~isfinite(D) || ~isfinite(n) || n<=0
    p = NaN; return; end
lambda = (sqrt(n) + 0.12 + 0.11/sqrt(n)) * D;
% Smirnov formula
terms = 1:100;
p = 2*sum((-1).^(terms-1).*exp(-2*(terms.^2)*(lambda^2)));
p = max(min(p,1),0);
end

