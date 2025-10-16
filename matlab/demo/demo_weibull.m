function demo_weibull()
%DEMO_WEIBULL Highlight Weibull fits for life/strength data.
setup_demo_paths();
rng(8);
scale = 520; shape = 2.1;
stressLife = wblrnd(scale, shape, 320, 1) + 40*rand(320,1);
stressLife = max(stressLife, 0);

opts = struct('orderBy','AD','computeP','none');
results = fitAllDistributions(stressLife, opts);

print_fit_table('Weibull Life Data Analysis', results, 5);
summaryTable = top_results_table(results, 5);
disp(summaryTable);

assignin('base','weibull_data', stressLife);
assignin('base','weibull_results', results);
assignin('base','weibull_summary', summaryTable);

bestModel = get_model_by_name(results(1).name);

figure('Name','Weibull Life Demo','NumberTitle','off');
subplot(1,2,1);
hist_pdf(stressLife, bestModel, results(1).theta, 30);
xlabel('Cycles to Failure'); ylabel('Density'); grid on;

subplot(1,2,2);
[F,t] = ecdf(stressLife);
stairs(t, F, 'b','LineWidth',1.2); hold on;
tt = linspace(min(stressLife), max(stressLife), 300);
Ft = bestModel.cdf(tt, results(1).theta);
plot(tt, Ft, 'r-','LineWidth',2);
hold off; grid on;
legend('Empirical CDF','Model CDF','Location','best');
xlabel('Cycles to Failure'); ylabel('Cumulative Probability');
title(sprintf('CDF Fit: %s', bestModel.name));
end

