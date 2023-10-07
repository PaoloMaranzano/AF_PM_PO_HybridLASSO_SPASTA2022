%% Check overlapping
function [Prop_fold] = check_prop_strf_cv(CV_part,Pollutant,Fold,Strat_var)

%%% Debugging
if 0
    Pollutant = 'PM10';
    Fold = 2;
    CV_part = cv_part;
    Strat_var = 'Tipology_rec';
end

poll = CV_part.poll;
p = find(ismember(poll, Pollutant));
if CV_part.cv_type == "StratBoot"
    % Stratified bootstrap case
    check_cv = CV_part.mat_fold{1,Fold}{1,p} == 1;
else
    % Other cases
    check_cv = CV_part.mat_fold{1,p} == Fold;
end
Prop_fold = cell2table(CV_part.tab{1,p});
Prop_fold.Properties.VariableNames = {'Strata','Obs_count','Share_pop'};
Prop_fold.(['Share_fold_' num2str(Fold)]) = zeros(CV_part.nStrata,1);
for Stratum = 1:CV_part.nStrata
    % Establish which stations belong to the stratum among those
    % which monitor the p-th pollutant (stations of interest)
    idx_st_poll = ismember(CV_part.ARPA_stats_reg.New_cod_stz,...
        CV_part.Coords_stats{1,p}.IDStat);
    idx_type_strat = ismember(CV_part.ARPA_stats_reg.(Strat_var),...
        CV_part.Strata(Stratum));
    SOI = find(idx_st_poll & idx_type_strat);
    idx_st_type_strat = idx_st_poll & idx_type_strat;
    idx_st_type_strat = idx_st_type_strat(idx_st_poll);
    % idx_st_type_strat = ismember(CV_part.Coords_stats{1,p}.IDStat,...
    %     categorical(string("stz_") + string(num2str(SOI,'%02d'))));
    % Share of stratum s-th in j-th fold
    sh_sample = sum(check_cv(idx_st_type_strat,:),'all') / sum(check_cv,'all');
    Prop_fold.(['Share_fold_' num2str(Fold)])(Stratum) = sh_sample*100;
end
return;

end
