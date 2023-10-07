%% CV partinioning: stratified K-fold cross-validation
function [cv_part] = CVpart_strat_Kfold_CondStat(DSTEM_str_ground,Strat_var,K)

%%% Debugging
if 0
    clear cv_part
    DSTEM_str_ground = Ground;
    Strat_var = 'Tipology_rec';
    K = 10;
end

% CV type
cv_part.cv_type = 'StratKFold';
% Folds
cv_part.K = K;
% Stations registry
cv_part.ARPA_stats_reg = DSTEM_str_ground.ARPA_stats_reg;
% Determine distinct strata
cv_part.Strata = unique(DSTEM_str_ground.ARPA_stats_reg.(Strat_var));
% Count distinct strata
cv_part.nStrata = size(cv_part.Strata,1);
% Pollutants
cv_part.poll = DSTEM_str_ground.poll;

%%% Generating partitions (K-folds)
for p = 1:length(cv_part.poll)
    %%% Store sites coordinates (necessary for check CV)
    cv_part.Coords_stats{1,p} = DSTEM_str_ground.(['Coords_' cv_part.poll{p}]);
    %%% Store useful quantities
    % Observed (original) response data
    cv_part.Y{1,p} = DSTEM_str_ground.([cv_part.poll{p}]);
    % Matrix with folder division
    cv_part.mat_fold{1,p} = DSTEM_str_ground.([cv_part.poll{p}]);
    % Non-missing observations index
    cv_part.not_nan_Y{1,p} = ~isnan(cv_part.Y{1,p});
    % Identify which stations monitor the p-th pollutant
    idx_st_poll = ismember(cv_part.ARPA_stats_reg.New_cod_stz,...
        DSTEM_str_ground.(['Coords_' cv_part.poll{p}]).IDStat);
    %%% Stratification conditioning to stations (rows)
    for row = 1:size(cv_part.mat_fold{1,p},1)
        % Non-missing values among the selected stations
        nonnan_y = ~isnan(cv_part.mat_fold{1,p}(row,:));
        % Count the non-missing values
        lng = length(find(nonnan_y));
        % Assign to all non-missing values in r-th row a value in 1 to K
        v = [repmat(1:K, 1, floor(lng/K)) , 1:(lng-floor(lng/K)*K)];
        % Random assignment to a folder (re-shuffle) by row
        cv_part.mat_fold{1,p}(row,nonnan_y) = datasample(v,lng,'Replace',false);
    end
    % Identify which stations monitor the p-th pollutant
    cv_part.idx_st_poll{1,p} = idx_st_poll;
    % Share of each stratum among the non-missing observed obs
    cv_part.tab{1,p} = repelem(cv_part.ARPA_stats_reg.(Strat_var)(cv_part.idx_st_poll{1,p},:),...
        sum(cv_part.not_nan_Y{1,p},2));
    cv_part.tab{1,p} = tabulate(string(cv_part.tab{1,p}));
    cv_part.tab{1,p}(find([cv_part.tab{1,p}{:,2}] == 0),:) = [];
end

%%% Splitting datasets into K-folds
for fold = 1:K
    for p = 1:length(cv_part.poll)
        cv_part.Y_cv{1,fold}{1,p} = cv_part.Y{1,p};
        % Identify candidate cells (non missing) for the p-th pollutant
        obs_fold = cv_part.mat_fold{1,p} == fold;
        % Store candidate cells coords (row and column)
        [r, c] = find(obs_fold);
        rc = [r,c];
        % Assign NaN to the stratified random sample
        for i = 1:size(rc,1)
            cv_part.Y_cv{1,fold}{1,p}(rc(i,1),rc(i,2)) = NaN;
        end
        % Store test set coords (row and col)
        cv_part.test_coords{1,fold}{1,p} = rc;
    end
end

return;

end