function demo_demand()
%DEMO_DEMAND Evaluate discrete demand models for operations planning.
setup_demo_paths();
rng(4);
demand = poissrnd(4.2, 260, 1) + binornd(2, 0.25, 260, 1);
demand = max(demand, 0);

opts = struct('orderBy','CS','computeP','none');
results = fitAllDistributions(demand, opts);

print_fit_table('Discrete Demand Counts', results, 5);
summaryTable = top_results_table(results, 5);
disp(summaryTable);

assignin('base','demand_counts', demand);
assignin('base','demand_results', results);
assignin('base','demand_summary', summaryTable);

bestModel = get_model_by_name(results(1).name);
kVals = (floor(min(demand)):ceil(max(demand))).';

figure('Name','Demand Model Fits','NumberTitle','off');
subplot(1,2,1);
histogram(demand, 'Normalization','probability','BinMethod','integers'); hold on;
pmf = bestModel.pdf(kVals, results(1).theta);
stem(kVals, pmf, 'r','LineWidth',1.5,'Marker','none');
hold off; grid on;
xlabel('Daily Demand'); ylabel('Probability');
legend('Empirical','Model PMF','Location','best');

subplot(1,2,2);
[F,x] = ecdf(demand);
stairs(x, F, 'b','LineWidth',1.2); hold on;
Fc = bestModel.cdf(kVals, results(1).theta);
stairs(kVals, Fc, 'r','LineWidth',1.5);
hold off; grid on;
xlabel('Daily Demand'); ylabel('Cumulative Probability');
legend('Empirical CDF','Model CDF','Location','best');
title(sprintf('Best Fit: %s', bestModel.name));
end

