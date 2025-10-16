function demo_gpd()
%DEMO_GPD Model threshold exceedances with the Generalized Pareto distribution.
setup_demo_paths();
rng(10);

threshold = 35;
gpdModel = dists.gpd();
exceedances = gpdModel.rnd(320, [5.2, 0.18]);
exceedances = exceedances + abs(0.4*randn(size(exceedances)));
tailSamples = threshold + exceedances;

opts = struct('orderBy','KS','computeP','none');
results = fitAllDistributions(exceedances, opts);

print_fit_table('Generalized Pareto Tail Modeling', results, 5);
summaryTable = top_results_table(results, 5);
disp(summaryTable);

assignin('base','gpd_threshold', threshold);
assignin('base','gpd_exceedances', exceedances);
assignin('base','gpd_tail_samples', tailSamples);
assignin('base','gpd_results', results);
assignin('base','gpd_summary', summaryTable);

bestModel = get_model_by_name(results(1).name);

figure('Name','GPD Threshold Exceedance Demo','NumberTitle','off');
subplot(1,2,1);
hist_pdf(exceedances, bestModel, results(1).theta, 25);
xlabel('Exceedance over Threshold'); ylabel('Density'); grid on;

tailProb = linspace(0, max(exceedances)*1.05, 200);
subplot(1,2,2);
S_emp = 1 - ecdf_step(exceedances, tailProb);
S_model = 1 - bestModel.cdf(tailProb, results(1).theta);
plot(tailProb, S_emp, 'b', 'LineWidth',1.5); hold on;
plot(tailProb, S_model, 'r', 'LineWidth',2);
hold off; grid on;
legend('Empirical Survival','Model Survival','Location','best');
xlabel('Exceedance'); ylabel('Tail Probability');
title(sprintf('Tail Survival: %s', bestModel.name));
end

function S = ecdf_step(x, grid)
[F, xi] = ecdf(x);
S = zeros(size(grid));
for i = 1:numel(grid)
    idx = find(xi <= grid(i), 1, 'last');
    if isempty(idx), S(i) = 0; else, S(i) = F(idx); end
end
end

