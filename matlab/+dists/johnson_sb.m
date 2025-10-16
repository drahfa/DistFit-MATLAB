function m = johnson_sb()
%JOHNSON_SB Placeholder for Johnson SB (bounded) distribution.
m.name = 'Johnson SB';
m.type = 'cont';
m.nParams = 4; % gamma, delta, xi, lambda
m.paramNames = {'gamma','delta','xi','lambda'};
m.support = [-inf, inf];
m.start = @(x) [0,1,min(x), max(x)-min(x)+eps];
m.fit = @(x,~) error('Johnson SB fitting not yet implemented');
m.pdf = @(x,th) NaN(size(x));
m.cdf = @(x,th) NaN(size(x));
m.rnd = [];
end
