function run_demo()
%RUN_DEMO Fit sample data files in ../data and print rankings.
root = fileparts(mfilename('fullpath'));
dataDir = fullfile(root, '..', 'data');
addpath(fullfile(root,'utils'));
addpath(fullfile(root,'gof'));
addpath(fullfile(root,'plot'));

files = { 'normal.txt', 'exp.txt', 'discrete.txt' };
for i = 1:numel(files)
    fn = fullfile(dataDir, files{i});
    if ~isfile(fn)
        fprintf('Missing file: %s\n', fn);
        continue;
    end
    x = read_numeric_file(fn);
    fprintf('\n=== %s ===\n', files{i});
    R = fitAllDistributions(x, struct('orderBy','KS'));
    print_results(R, 5);
end
end

function print_results(R, k)
if nargin<2, k = numel(R); end
k = min(k, numel(R));
fprintf('%-12s  %-20s  %-8s  %-8s  %-8s\n','Rank','Name','KS','AD','ChiSq');
for i=1:k
    r = R(i);
    fprintf('#%-11d  %-20s  %-8.4g  %-8.4g  %-8.4g\n', i, r.name, r.KS, r.AD, r.CS);
end
end
