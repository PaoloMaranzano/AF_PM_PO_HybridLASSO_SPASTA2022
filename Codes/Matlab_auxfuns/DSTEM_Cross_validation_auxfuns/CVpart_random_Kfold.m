%% CV partinioning: random K-fold cross-validation
function [cv_part] = CVpart_random_Kfold(DSTEM_str_ground,K)

%%% Debugging
if 0
    clear cv_part
    DSTEM_str_ground = Ground;
    K = 10;
end

% CV type
cv_part.cv_type = 'RandomKFold';
% Folds
cv_part.K = K;
% Stations registry
cv_part.ARPA_stats_reg = DSTEM_str_ground.ARPA_stats_reg;
% Response variables name
cv_part.poll = DSTEM_str_ground.Y_names;

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
    % Count the non-missing values
    lng = length(find(nonnan_y));
    % Assign to all non-missing values in r-th row a value in 1 to K
    v = [repmat(1:K, 1, floor(lng/K)) , 1:(lng-floor(lng/K)*K)];
    % Random assignment to a folder (re-shuffle) by row
    cv_part.mat_fold{1,p}(nonnan_y) = datasample(v,lng,'Replace',false);
    % Identify which stations monitor the p-th pollutant
    cv_part.idx_st_poll{1,p} = idx_st_poll;
end

%%% guardare riga 92 extract
%%% y_mat è trasposta --> cv_part.Y = y_mat va ritraposta prima di (:)

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
        % Store Y for functional case (long for data_table)
		if ismember(DSTEM_str_ground.model_type,"f-HDGM")
            % mat = st x T --> long = mat'(:)
            Y_cv_t = cv_part.Y_cv{1,fold}{1,p}'; % 24/03 ' sbagliato
            Y_long_cv = Y_cv_t(:);
            cv_part.Y_long_cv{1,fold}{1,p} = Y_long_cv;
            Y_split = splitapply( @(x){x}, Y_long_cv, DSTEM_str_ground.data_long.gId );
            % Y_split = splitapply( @(x){x}, Y_long_cv, gId );
            Y_fun_cv = cellfun(@transpose,Y_split,'UniformOutput',false);
            cv_part.Y_fun_cv{1,fold}{1,p} = Y_fun_cv;
        end
    end
end

% Store Y for functional case (long for data_table)
if ismember(DSTEM_str_ground.model_type,"f-HDGM")
    for p = 1:length(cv_part.poll)
        Y_t = cv_part.Y{1,p}';
        Y_long = Y_t(:);
        cv_part.Y_long{1,p} = Y_long;
        Y_split = splitapply( @(x){x}, Y_long, DSTEM_str_ground.data_long.gId );
        Y_fun = cellfun(@transpose,Y_split,'UniformOutput',false);
        cv_part.Y_fun{1,p} = Y_fun;
    end
end

return;

end