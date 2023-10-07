%% CV partinioning: stratified K-fold cross-validation
function cv_part = CVpart_strat_Kfold(DSTEM_str_ground,Strat_var,K)

%%% Debugging
if 0
    clear cv_part
    DSTEM_str_ground = Ground;
    Strat_var = 'Tipology_rec';
    K = 10;
end

% CV type
cv_part.cv_type = 'StratKFold';
% Runs
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
    %%% Stratification
    for Stratum = 1:cv_part.nStrata
        % Establish which stations belong to the stratum among those which monitor 
        % the p-th pollutant (stations of interest)
        idx_type_strat = ismember(cv_part.ARPA_stats_reg.(Strat_var),...
            cv_part.Strata(Stratum));
        SOI = find(idx_st_poll & idx_type_strat);
        % idx_st_type_strat = ismember(cv_part.Coords_stats{1,p}.IDStat,...
        %      categorical(string("stz_") + string(num2str(SOI,'%02d'))));
        idx_st_type_strat = idx_st_poll & idx_type_strat;
        idx_st_type_strat = idx_st_type_strat(idx_st_poll);
        % Support matrix which shows selected stations
        mat_st = repmat(idx_st_type_strat,1,size(cv_part.Y{1,p},2));
        % Non-missing values among the selected stations
        nonnan_y = ~isnan(cv_part.mat_fold{1,p});
        % Count the non-missing values
        lng = length(find(nonnan_y(idx_st_type_strat,:)));
        % Assign to all the non-missing values a value in 1 to K
        v = [repmat(1:K, 1, floor(lng/K)) , 1:(lng-floor(lng/K)*K)];
        % Random assignment to a folder (re-shuffle)
        cv_part.mat_fold{1,p}(find(nonnan_y .* mat_st)) = datasample(v,lng,'Replace',false);
    end
    % Identify which stations monitor the p-th pollutant
    cv_part.idx_st_poll{1,p} = idx_st_poll;
    % Share of each stratum among the non-missing observed obs
    cv_part.tab{1,p} = repelem(cv_part.ARPA_stats_reg.(Strat_var)(cv_part.idx_st_poll{1,p},:),...
        sum(cv_part.not_nan_Y{1,p},2));
    % cv_part.tab{1,p} = tabulate(string(cv_part.tab{1,p}));
    cv_part.tab{1,p} = tabulate(cv_part.tab{1,p});
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