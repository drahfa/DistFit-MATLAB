function demo_exploratory()
%DEMO_EXPLORATORY Rank candidate models for a mixed continuous sample.
setup_demo_paths();
rng(1); % reproducible sample
x = [1.2*randn(320,1) + 1.0; 0.6*randn(180,1) - 0.8];

opts = struct('orderBy','KS','computeP','none');
results = fitAllDistributions(x, opts);

print_fit_table('Exploratory Data Analysis', results, 5);
summaryTable = top_results_table(results, 5);
disp(summaryTable);

assignin('base','eda_data',x);
assignin('base','eda_results',results);
assignin('base','eda_summary',summaryTable);

bestModel = get_model_by_name(results(1).name);

figure('Name','EDA Distribution Fits','NumberTitle','off');
subplot(1,2,1);
hist_pdf(x, bestModel, results(1).theta);
xlabel('Observation'); ylabel('Density'); grid on;

subplot(1,2,2);
ecdf_cdf(x, bestModel, results(1).theta);
xlabel('Observation'); ylabel('Cumulative Probability'); grid on;
end

