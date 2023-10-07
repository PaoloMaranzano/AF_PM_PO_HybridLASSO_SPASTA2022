%% %%%%% Stratified K-fold cross-validation for spatio-temporal models
% Testing observations are randomly sampled and assigned to a fold aiming
% at keeping the proportions of a categorical variable. Each fold has to
% represent the composition of the original dataset. Each fold is composed
% inter_lock_type = 1;
% run('..\..\..\Data\DSTEM_data\Reshape_to_HDGM.m');
clearvars -except crossval_step data Ground log_transform poll standardization
debug_cv = 1;
save_cv = 0;
% poll_list = {'NO2'};
% p = 1;
poll = {'NO2'};
standardization = 1;
log_transform = 1;

k = 2;
MSE = cell(length(poll),k);
RMSE = cell(length(poll),k);
MAE = cell(length(poll),k);
R2 = cell(length(poll),k);

%%%%% Objective: draw a sample from each row (station)
%%% For each row (station):
% Non-missing observations are assigned randomly to a fold by assigning
% them a value between 1 and k.
% If T(nonnan)/k is even in each row, all the folds are balanced: same
% number of non-missing obs.
% If T(nonnan)/k is odd in each row, folds are unbalanced: different 
% number of non-missing obs.
%%% !!! For each row, the number of non-missing value can be different,
%%% hence the number of observations in each folder can differ. In
%%% particular, at increasing k corresponds a decreasing number of obs.
%%% The variation is around 1%
for p = 1:length(poll)
    cv_part.Y{1,p} = Ground.([poll{p}]);
    cv_part.mat_fold{1,p} = Ground.([poll{p}]);
    for row = 1:size(cv_part.mat_fold{1,p},1)
        nonnan_y = ~isnan(cv_part.mat_fold{1,p}(row,:));
        l = length(find(nonnan_y));
        v = [repmat(1:k, 1, floor(l/k)) , 1:(l-floor(l/k)*k)];
        cv_part.mat_fold{1,p}(row,nonnan_y) = datasample(v,l,'Replace',false);
    end
end

% Determine distinct strata
DistinctStrata = unique(Ground.ARPA_stats_reg.Tipology_rec);
% Count distinct strata
nDistinctStrata = size(DistinctStrata,1);

parfor run = 1:k
    %% Computation begin
    date_begin = datetime('now')
    %% Loading data
    % load('..\..\Data\DSTEM_data\HDGM_data_SPASTA.mat')
    inter_lock_type = Ground.InterLockType;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%      Data  building     %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    n = cell(length(poll),1);
    obj_stem_gridlist_p = stem_gridlist();
    ground = struct();
    for p = 1:length(poll)
        %%% Dependent variables
        ground.Y{p} = Ground.([poll{p}]);
        ground.Y_name{p} = [poll{p} ' ground'];
        n{p,1} = size(ground.Y{p}, 1);
        %%% Loading covariates for the selected pollutants at each monitoring stations
        ground.X_beta{p} = Ground.(['X_' poll{p}]);
        ground.X_beta_name{p} = Ground.vars_names;
        ground.X_z{p} = ones(n{p}, 1);
        ground.X_z_name{p} = {['constant_' poll{p}]};
        %%% Coordinates grid
        % stem_grid: contains all the information related to the sampling
        % locations of a single variable
        ground.coordinates{p} = [Ground.(['Coords_' poll{p}]).Latitude, ...
            Ground.(['Coords_' poll{p}]).Longitude];
        obj_stem_grid = cell(length(poll),1);
        obj_stem_grid{p} = stem_grid(ground.coordinates{p}, 'deg', 'sparse', 'point');
        % stem_gridlist: collector of the stem_grid objects for all the variables;
        obj_stem_gridlist_p.add(obj_stem_grid{p});
    end
    T = size(ground.Y{1}, 2);
    
    for p = 1:length(poll)
        ground.Y_cv{1,p} = ground.Y{1,p};
        % Identify non-missing values for the p-th pollutant
        not_nan_poll = ~isnan(ground.Y_cv{1,p});
        % Count of non-missing values for the p-th pollutant
        cnt_not_nan_poll = sum(not_nan_poll,'all');
        % Test set size
        test_size = floor(cnt_not_nan_poll/k);
        % Identify which stations monitor the p-th pollutant
        idx_st_poll = ismember(Ground.ARPA_stats_reg.New_cod_stz,Ground.(['Coords_' poll{p}]).IDStat);
        % Generate partitions (weights) by stations type by computing the
        % sharing of each strata among the non-missing observed obs
        tab = repelem(Ground.ARPA_stats_reg.Tipology_rec(idx_st_poll,:), ...
            sum(not_nan_poll,2));
        tab = tabulate(tab);
        % Identify candidate cells (non missing) for the p-th pollutant
        mat_rnd_not_nan = cv_part.m{1,p} .* not_nan_poll;
        not_nan_poll = mat_rnd_not_nan == run;
        % Store candidate cells coords (row and column)
        [row_not_nan, col_not_nan] = find(not_nan_poll);
        rc = [row_not_nan, col_not_nan];
        % Stratification
        for Stratum = 1:nDistinctStrata
            % Establish which stations belong to the stratum among those
            % which monitor the p-th pollutant (stations of interest)
            idx_type_strat = ismember(Ground.ARPA_stats_reg.Tipology_rec,DistinctStrata(Stratum));
            SOI = find(idx_st_poll & idx_type_strat);
            % Count the number of station in the stratum
            nSOI = length(SOI);
            % Sample the stations to use as k-fold
                % Set r = 1 to draw randomly from all the stations 
                % (stratified sampling)
                % Set 0<r<1 to draw randomly from a random subset of the 
                % stations (cluster sampling)
            r = 1;
            idx_st_type_strat = ismember(Ground.(['Coords_' poll{p}]).IDStat,...
                categorical(string("stz_") + string(num2str(SOI,'%02d'))));
            idx_S_sel = find(idx_st_type_strat);
            S_sel = datasample(idx_S_sel,round(r*nSOI),'Replace',false);
            % Identify candidate cells (non missing) in the stratum
            st = rc(ismember(rc(:,1),S_sel),:);
            % Stratified test size (weight of the stratum)
            perc = tab{tab(:,1) == DistinctStrata(Stratum),3}/100;
            % Random sample from the candidates in the stratum
            rc_sel = datasample(st,min(floor(test_size*perc),size(st,1)),...
                'Replace',false);
            % Assign NaN to the stratified random sample
            for i = 1:size(rc_sel,1)
                ground.Y_cv{1,p}(rc_sel(i,1),rc_sel(i,2)) = NaN;
            end
        end
        % Store test set coords (row and col) and strata weights
        [r,c] = find(isnan(ground.Y_cv{1,p}) - isnan(ground.Y{1,p}));
        ground.test_coords.(poll{p}) = table(r);
        ground.test_coords.(poll{p}).c = c;
        ground.tab.(poll{p}) = tab;
    end
    % Store Y original data
    if run == 1
        y_original{1,run} = ground.Y;
    end
    % Store Y training set
    y_train_cv{1,run} = ground.Y_cv;
    % Store Y training set coordinates
    coord_y_test{1,run} = ground.test_coords;
    % Store Y training tabulate
    tab_y{1,run} = ground.tab;
    
    if debug_cv == 0
        %% %%%%% STEM settings
        % stem_varset: contains the observed data of all the variables and the
        % loading coefficients;
        Y_obs = struct();
        Y_obs = ground.Y;
        test_coords = ground.test_coords;
        obj_stem_varset_p = stem_varset(ground.Y_cv, ground.Y_name, [], [], ...
            ground.X_beta, ground.X_beta_name, ...
            ground.X_z, ground.X_z_name);
        
        % datestamp: contains the information related to the date and time of the
        % observations;
        obj_stem_datestamp = stem_datestamp('01-01-2017 00:00','31-10-2020 00:00',T);
        
        % HDGM STEM_data object
        % Shapefile
        shape = [];
        obj_stem_modeltype = stem_modeltype('HDGM');
        obj_stem_data = stem_data(obj_stem_varset_p, obj_stem_gridlist_p, ...
            [], [], obj_stem_datestamp, [], ...
            obj_stem_modeltype, shape);
        
        
        %% %%%%% DSTEM_par
        % DSTEM_par object creation: contains the structure and the values of the
        % model parameters;
        obj_stem_par_constraints=stem_par_constraints();
        % time_diagonal: used to specify if the matrices G and Sigma_eta are
        % diagonal or not
        obj_stem_par_constraints.time_diagonal=0;
        obj_stem_par = stem_par(obj_stem_data, 'exponential',obj_stem_par_constraints);
        % DSTEM_model object creation
        obj_stem_model = stem_model(obj_stem_data, obj_stem_par);
        % clear ground
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%      Data transform      %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if log_transform == 1
            obj_stem_model.stem_data.log_transform;
        end
        if standardization == 1
            obj_stem_model.stem_data.standardize;
        end
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%      Model estimation      %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Starting values
        obj_stem_par.beta = obj_stem_model.get_beta0();
        obj_stem_par.theta_z = 100;         % Kilometers
        obj_stem_par.v_z = [
            1 0.6 0.6;
            0.6 1 0.9
            0.6 0.9 1];  % Cross-correlation between multiple variables
        obj_stem_par.sigma_eta = diag([0.2 0.2 0.2]);
        obj_stem_par.G = diag([0.8 0.8 0.8]);
        obj_stem_par.sigma_eps = diag([0.3 0.3 0.3]);
        obj_stem_model.set_initial_values(obj_stem_par);
        
        % Parameters estimation
        obj_stem_EM_options = stem_EM_options();
        obj_stem_EM_options.exit_tol_loglike = 0.0005;
        obj_stem_EM_options.exit_tol_par = 0.0005;
        obj_stem_EM_options.max_iterations = 300;
        date_estimation_begin = datetime('now')
        obj_stem_model.EM_estimate(obj_stem_EM_options);
        date_estimation_end = datetime('now')
        
        y_oos = struct();
        for p = 1:length(poll)
            test_coords.(poll{p}) = table2array(test_coords.(poll{p}));
            for i = 1:size(test_coords.(poll{p}),1)
                y_oos.(poll{p})(i,1) = test_coords.(poll{p})(i,1);
                y_oos.(poll{p})(i,2) = test_coords.(poll{p})(i,2);
                y_oos.(poll{p})(i,3) = Y_obs{1,1}(test_coords.(poll{p})(i,1),test_coords.(poll{p})(i,2));
                y_oos.(poll{p})(i,4) = obj_stem_model.stem_EM_result.y_hat{1,1}(test_coords.(poll{p})(i,1),...
                    test_coords.(poll{p})(i,2));
            end
            y_oos.(poll{p}) = array2table(y_oos.(poll{p}));
            y_oos.(poll{p}).Properties.VariableNames = {'row','col','y_obs','y_hat'};
            y_oos.MSE{p,1} = nanmean((y_oos.(poll{p}).y_hat - y_oos.(poll{p}).y_obs).^2);
            y_oos.RMSE{p,1} = sqrt(y_oos.MSE{p,1});
            y_oos.MAE{p,1} = nanmean(abs(y_oos.(poll{p}).y_hat - y_oos.(poll{p}).y_obs));
            y_oos.R2{p,1} = 1 - y_oos.MSE{p,1} ./ nanvar(y_oos.(poll{p}).y_obs);
        end
        
        MSE(:,run) = y_oos.MSE;
        RMSE(:,run) = y_oos.RMSE;
        MAE(:,run) = y_oos.MAE;
        R2(:,run) = y_oos.R2;
    end
end

cv_output.MSE = MSE;
cv_output.RMSE = RMSE;
cv_output.MAE = MAE;
cv_output.R2 = R2;

if save_cv == 1
    save('strat_kfold_cv.mat');
end











