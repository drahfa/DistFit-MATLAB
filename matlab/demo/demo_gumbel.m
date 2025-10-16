function demo_gumbel()
%DEMO_GUMBEL Explore extreme-value fits for block maxima data.
setup_demo_paths();
rng(7);
mu = 38; sigma = 5.5;
base = dists.gumbel().rnd(360, [mu, sigma]);
noise = 0.6*randn(size(base));
blockMax = base + noise;

opts = struct('orderBy','KS','computeP','none');
results = fitAllDistributions(blockMax, opts);

print_fit_table('Gumbel Extreme Value Analysis', results, 5);
summaryTable = top_results_table(results, 5);
disp(summaryTable);

assignin('base','gumbel_data', blockMax);
assignin('base','gumbel_results', results);
assignin('base','gumbel_summary', summaryTable);

bestModel = get_model_by_name(results(1).name);

figure('Name','Gumbel Extreme Value Demo','NumberTitle','off');
subplot(1,2,1);
hist_pdf(blockMax, bestModel, results(1).theta, 25);
xlabel('Annual Maximum (mm)'); ylabel('Density'); grid on;

subplot(1,2,2);
qqplot_model(blockMax, bestModel, results(1).theta);
grid on;
title(sprintf('Q-Q Plot against %s', bestModel.name));
end

