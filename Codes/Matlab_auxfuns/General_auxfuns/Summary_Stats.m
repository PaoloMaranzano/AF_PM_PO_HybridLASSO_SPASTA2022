function [stats_tab] = Summary_Stats(var)

namestr = inputname(1);
namestr = string(namestr);

n = length(var);
minn = min(var);
maxx = max(var);
m = nanmean(var);
s2 = nanvar(var);
s = sqrt(s2);
sem = s/sqrt(n);
VC = s/m;
sk = skewness(var);
k = kurtosis(var);
q1 = quantile(var,0.01);
q5 = quantile(var,0.05);
q25 = quantile(var,0.25);
q50 = quantile(var,0.50);
q75 = quantile(var,0.75);
q95 = quantile(var,0.95);
q99 = quantile(var,0.99);

stats_tab = [n,minn,q1,q5,q25,m,q50,q75,q95,q99,maxx,s2,s,sem,VC,sk,k];
stats_tab = array2table(stats_tab);
stats_tab.Properties.RowNames = namestr;
stats_tab.Properties.VariableNames = {'n','min','q1','q5','q25','mean','median',...
    'q75','q95','q99','max','s2','s','sem','VC','skew','kurt'};


end