function [ActiveSets_opt,red_ASi,avg_dist_all,dist_m_full] = DSTEM_CV_AS_reduction(...
    ActiveSets,cv_metrics,metric)

if 0
    cv_metrics = randK_CV_metrics;
    metric = 'RMSE';
    ActiveSets = ActiveSets;
end

% lockdown and restart in AS
ActiveSets_lockrest = ActiveSets(contains(ActiveSets.Variable,{'lock1','restart1'}),:);
% AS which contains at least one of lock or restart for each pollutant
ActiveSets_lockrest_cnt = grpstats(ActiveSets_lockrest(:,2:end),'Response_idx','sum');

% Useful quantities
posAS = find(contains(ActiveSets.Properties.VariableNames,'Active'));
nAS = sum(contains(ActiveSets.Properties.VariableNames,'Active'));
nY = max(ActiveSets.Response_idx);

for ASi = 1:nAS
    incl_all_p(ASi) = sum(table2array(ActiveSets_lockrest_cnt(:,posAS(ASi))) >= [1 1 1]') == nY;    
end
red_ASi = find(incl_all_p);

for p = 1:nY
    metr_int = cell2mat(cv_metrics.([metric '_ASi']){1,p});
    dist_m_full{p} = metr_int - metr_int(:,nAS);
    avg_dist(p,:) = mean(dist_m_full{p},1);
end
avg_dist_all = avg_dist;
avg_dist = avg_dist(:,red_ASi);
for p = 1:nY
    pos_p{p} = find(avg_dist(p,:) <= 0);
    opt_AS_p{p} = pos_p{p}(avg_dist(pos_p{p}) == min(avg_dist(pos_p{p}))) + 2;
end
avg_dist_gen = mean(avg_dist,1);
pos_avg = find(avg_dist_gen <= 0);
opt_avg = pos_avg(avg_dist_gen(pos_avg) == min(avg_dist_gen(pos_avg))) + 2;

ActiveSets_red = ActiveSets(:,[1 , 2 , red_ASi]);
ActiveSets_opt = ActiveSets_red(:,[1 , 2, [opt_AS_p{:}], opt_avg]);




return;
end