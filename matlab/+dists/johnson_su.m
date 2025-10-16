function m = johnson_su()
%JOHNSON_SU Placeholder for Johnson SU (unbounded) distribution.
m.name = 'Johnson SU';
m.type = 'cont';
m.nParams = 4; % gamma, delta, xi, lambda
m.paramNames = {'gamma','delta','xi','lambda'};
m.support = [-inf, inf];
m.start = @(x) [0,1,mean(x,'omitnan'), std(x,0,'omitnan')+eps];
m.fit = @(x,~) error('Johnson SU fitting not yet implemented');
m.pdf = @(x,th) NaN(size(x));
m.cdf = @(x,th) NaN(size(x));
m.rnd = [];
end
