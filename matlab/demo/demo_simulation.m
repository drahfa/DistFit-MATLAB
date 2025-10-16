function demo_simulation()
%DEMO_SIMULATION Fit a service-time model and generate synthetic samples.
setup_demo_paths();
rng(5);
serviceTimes = lognrnd(1.9, 0.38, 360, 1) + 3;

opts = struct('orderBy','KS','computeP','none');
results = fitAllDistributions(serviceTimes, opts);

print_fit_table('Simulation via Best-Fit Model', results, 5);
summaryTable = top_results_table(results, 5);
disp(summaryTable);

assignin('base','simulation_data', serviceTimes);
assignin('base','simulation_results', results);
assignin('base','simulation_summary', summaryTable);

bestModel = get_model_by_name(results(1).name);

if isfield(bestModel, 'rnd') && ~isempty(bestModel.rnd)
    syntheticSamples = bestModel.rnd(1000, results(1).theta);
else
    warning('Selected model %s does not provide an RNG; falling back to empirical bootstrap.', bestModel.name);
    idx = randi(numel(serviceTimes), 1000, 1);
    syntheticSamples = serviceTimes(idx);
end

metrics = {@mean, @std, @(y) prctile(y, 95)};
metricNames = {'Mean','StdDev','P95'};
observedStats = cellfun(@(f) f(serviceTimes), metrics);
syntheticStats = cellfun(@(f) f(syntheticSamples), metrics);
comparisonTable = table(metricNames.', observedStats.', syntheticStats.', ...
    'VariableNames', {'Metric','Observed','Synthetic'});
disp(comparisonTable);

assignin('base','synthetic_samples', syntheticSamples);
assignin('base','simulation_metrics', comparisonTable);

figure('Name','Simulation Fit and Synthetic Data','NumberTitle','off');
subplot(1,2,1);
hist_pdf(serviceTimes, bestModel, results(1).theta, 25);
xlabel('Service Time'); ylabel('Density'); grid on;

subplot(1,2,2);
histogram(serviceTimes, 'Normalization','pdf', 'FaceAlpha',0.6); hold on;
histogram(syntheticSamples, 'Normalization','pdf', 'FaceAlpha',0.6);
hold off; grid on;
legend('Observed','Synthetic', 'Location','best');
xlabel('Service Time'); ylabel('Density');
title(sprintf('Synthetic Sampling via %s', bestModel.name));
end
