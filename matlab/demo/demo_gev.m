function demo_gev()
%DEMO_GEV Fit generalized extreme value models to block maxima data.
setup_demo_paths();
rng(9);
mu = 55; sigma = 8; kappa = 0.18;
gevModel = dists.gev();
maxima = gevModel.rnd(400, [mu, sigma, kappa]);
maxima = maxima + 0.8*randn(size(maxima));

opts = struct('orderBy','KS','computeP','none');
results = fitAllDistributions(maxima, opts);

print_fit_table('Generalized Extreme Value Analysis', results, 5);
summaryTable = top_results_table(results, 5);
disp(summaryTable);

assignin('base','gev_data', maxima);
assignin('base','gev_results', results);
assignin('base','gev_summary', summaryTable);

bestModel = get_model_by_name(results(1).name);

figure('Name','GEV Extreme Value Demo','NumberTitle','off');
subplot(1,2,1);
hist_pdf(maxima, bestModel, results(1).theta, 30);
xlabel('Block Maximum'); ylabel('Density'); grid on;

subplot(1,2,2);
ecdf_cdf(maxima, bestModel, results(1).theta);
xlabel('Block Maximum'); ylabel('Cumulative Probability'); grid on;
title(sprintf('CDF Overlay: %s', bestModel.name));
end
