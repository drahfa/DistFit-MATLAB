function demo_reliability()
%DEMO_RELIABILITY Compare lifetime models for accelerated failure data.
setup_demo_paths();
rng(2);
lifetimes = wblrnd(850, 2.4, 260, 1) + exprnd(120, 260, 1);
lifetimes = lifetimes + 150; % enforce realistic lower bound

opts = struct('orderBy','AD','computeP','none');
results = fitAllDistributions(lifetimes, opts);

print_fit_table('Reliability / Lifetime Data', results, 5);
summaryTable = top_results_table(results, 5);
disp(summaryTable);

assignin('base','reliability_data', lifetimes);
assignin('base','reliability_results', results);
assignin('base','reliability_summary', summaryTable);

bestModel = get_model_by_name(results(1).name);

figure('Name','Reliability Fits','NumberTitle','off');
subplot(1,2,1);
hist_pdf(lifetimes, bestModel, results(1).theta, 25);
xlabel('Time to Failure (hrs)'); ylabel('Density'); grid on;

subplot(1,2,2);
[F,t] = ecdf(lifetimes);
stairs(t, 1-F, 'b','LineWidth',1.2); hold on;
tt = linspace(min(lifetimes), max(lifetimes), 300);
surv = 1 - bestModel.cdf(tt, results(1).theta);
plot(tt, surv, 'r-','LineWidth',2);
hold off; grid on;
ylabel('Reliability'); xlabel('Time to Failure (hrs)');
legend('Empirical Survival','Model Survival','Location','best');
title(sprintf('Reliability Curve: %s', bestModel.name));
end

