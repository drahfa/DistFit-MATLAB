function ecdf_cdf(x, M, theta)
%ECDF_CDF Plot empirical CDF and overlay model CDF
x = x(:); x = x(isfinite(x));
[f,xx] = ecdf(x);
stairs(xx, f, 'b'); hold on;
xx2 = linspace(min(x), max(x), 300);
F = M.cdf(xx2, theta);
plot(xx2, F, 'r-', 'LineWidth',2);
title(sprintf('%s CDF fit', M.name));
legend('Empirical','Model','Location','best');
hold off;
end

