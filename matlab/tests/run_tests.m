function run_tests()
%RUN_TESTS Validate core fitting logic on synthetic and sample data.
% Usage: run in MATLAB: run('matlab/tests/run_tests.m')

this = mfilename('fullpath');
testsDir = fileparts(this);
root = fileparts(testsDir); % matlab/
addpath(root);
addpath(fullfile(root,'utils'));

fprintf('Running tests in %s\n', testsDir);

%% Synthetic: Normal
rng(1);
mu = 30; sigma = 5;
x = mu + sigma*randn(10000,1);
R = fitAllDistributions(x, struct('orderBy','KS'));
assert_top(R, 'Normal');
th = R(strcmp({R.name},'Normal')).theta;
assert(abs(th(1)-mu) < 0.2, 'Normal mu off: got %.3f', th(1));
assert(abs(th(2)-sigma) < 0.2, 'Normal sigma off: got %.3f', th(2));
fprintf('Synthetic Normal: PASSED\n');

%% Synthetic: Exponential
rng(2);
muE = 12;
% exprnd fallback: inverse CDF
x = -muE*log(max(1e-12, 1-rand(10000,1)));
R = fitAllDistributions(x, struct('orderBy','KS'));
[idx, model] = find_model(R, 'Exponential');
assert(~isempty(idx), 'Exponential not in results');
assert(idx <= 2 || model.KS <= R(1).KS + 0.02, 'Exponential not near top (idx=%d)', idx);
th = R(strcmp({R.name},'Exponential')).theta;
assert(abs(th(1)-muE) < 0.5, 'Exponential mu off: got %.3f', th(1));
fprintf('Synthetic Exponential: PASSED\n');

%% Synthetic: Poisson (discrete)
rng(3);
lam = 5;
x = poiss_knuth(lam, 2000);
R = fitAllDistributions(x, struct('orderBy','KS'));
[idxP, modelP] = find_model(R, 'Poisson');
assert(~isempty(idxP), 'Poisson not in results');
assert(idxP <= 2 || modelP.KS <= R(1).KS + 0.03, 'Poisson not near top (idx=%d)', idxP);
th = R(strcmp({R.name},'Poisson')).theta;
assert(abs(th(1)-lam) < 0.5, 'Poisson lambda off: got %.3f', th(1));
fprintf('Synthetic Poisson: PASSED\n');

%% Provided data files
dataDir = fullfile(root,'..','data');
% normal.txt
xn = read_numeric_file(fullfile(dataDir,'normal.txt'));
Rn = fitAllDistributions(xn, struct('orderBy','KS'));
assert_top(Rn, 'Normal');
fprintf('Data normal.txt: PASSED\n');
% exp.txt (decimal comma)
xe = read_numeric_file(fullfile(dataDir,'exp.txt'));
Re = fitAllDistributions(xe, struct('orderBy','KS'));
[idxE, modelE] = find_model(Re, 'Exponential');
assert(~isempty(idxE), 'Exponential not in results for exp.txt');
assert(idxE <= 2 || modelE.KS <= Re(1).KS + 0.03, 'Exponential not near top for exp.txt (idx=%d)', idxE);
fprintf('Data exp.txt: PASSED\n');
% discrete.txt
xd = read_numeric_file(fullfile(dataDir,'discrete.txt'));
Rd = fitAllDistributions(xd, struct('orderBy','KS'));
[idxPD, modelPD] = find_model(Rd, 'Poisson');
assert(~isempty(idxPD), 'Poisson not in results for discrete.txt');
assert(idxPD <= 2 || modelPD.KS <= Rd(1).KS + 0.03, 'Poisson not near top for discrete.txt (idx=%d)', idxPD);
fprintf('Data discrete.txt: PASSED\n');

fprintf('\nAll tests PASSED.\n');

end

function assert_top(R, expected)
assert(~isempty(R), 'Empty results');
top = R(1).name;
assert(strcmpi(top, expected), 'Top model=%s, expected=%s', top, expected);
end

function [idx, model] = find_model(R, name)
idx = find(strcmpi({R.name}, name), 1, 'first');
if isempty(idx)
    model = struct();
else
    model = R(idx);
end
end

function x = poiss_knuth(lambda, n)
% Generate Poisson via Knuth algorithm (for moderate lambda)
if nargin<2, n = 1; end
x = zeros(n,1);
L = exp(-lambda);
for i=1:n
    k = 0; p = 1;
    while p > L
        k = k + 1;
        p = p * rand();
    end
    x(i) = k - 1;
end
end
