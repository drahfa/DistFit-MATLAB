function demo_finance()
%DEMO_FINANCE Fit heavy-tailed return distributions for risk analytics.
setup_demo_paths();
rng(3);
n = 620;
returns = 0.0004 + 0.0105*randn(n,1);
tailFlag = rand(n,1) < 0.08;
returns(tailFlag) = returns(tailFlag) + 0.035*trnd(4, sum(tailFlag), 1);

opts = struct('orderBy','KS','computeP','none');
results = fitAllDistributions(returns, opts);

print_fit_table('Financial Returns (Heavy Tails)', results, 5);
summaryTable = top_results_table(results, 5);
disp(summaryTable);

assignin('base','finance_returns', returns);
assignin('base','finance_results', results);
assignin('base','finance_summary', summaryTable);

bestModel = get_model_by_name(results(1).name);

figure('Name','Financial Return Fits','NumberTitle','off');
subplot(1,2,1);
hist_pdf(returns, bestModel, results(1).theta, 30);
xlabel('Return'); ylabel('Density'); grid on;

subplot(1,2,2);
qqplot_model(returns, bestModel, results(1).theta);
grid on;
end

