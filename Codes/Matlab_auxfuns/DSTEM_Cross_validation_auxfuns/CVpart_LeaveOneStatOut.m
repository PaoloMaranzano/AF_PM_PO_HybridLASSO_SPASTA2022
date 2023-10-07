%% CV partinioning: leave-one-station-out cross-validation
function [cv_part] = CVpart_LeaveOneStatOut(DSTEM_str_ground)

%%% Debugging
if 0
    clear cv_part
    DSTEM_str_ground = Ground;
end

% CV type
cv_part.cv_type = 'LeaveOneStatOut';
% Stations registry
cv_part.ARPA_stats_reg = DSTEM_str_ground.ARPA_stats_reg;
% Dependent variables
cv_part.poll = DSTEM_str_ground.poll;
% Folds
cv_part.K = size(cv_part.ARPA_stats_reg,1);

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
    % Non-missing values among the selected stations
    nonnan_y = cv_part.not_nan_Y{1,p};
    % Assign to each row (station) the corresponding code (1 to K)
    stat_lng = sum(idx_st_poll);
    stat_id = find(idx_st_poll);
    for s = 1:stat_lng
        cv_part.mat_fold{1,p}(s,nonnan_y(s,:)) = stat_id(s);
    end
    % Identify which stations monitor the p-th pollutant
    cv_part.idx_st_poll{1,p} = idx_st_poll;
end

%%% Splitting datasets into K-folds
for fold = 1:cv_part.K
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