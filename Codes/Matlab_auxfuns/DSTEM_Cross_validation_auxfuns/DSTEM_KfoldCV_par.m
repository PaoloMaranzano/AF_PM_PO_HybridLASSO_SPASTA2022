%% %%%%%%%%%% DSTEM_KfoldCV_par
%%% Evaluate K-fold CV performances of a DSTEM model using a parallelized internal for loop.
%%% See CVpart for avaiable CV schemes

function [CV_perf_metrics] = DSTEM_KfoldCV_par(cv_part,DSTEM_model,ActiveSets,ASi_idx)

if 0
    DSTEM_model = obj_stem_model_full;
    ASi_idx = 18;
    cv_part = cv_part_LOSO;
end
ASi_cols = find(contains(ActiveSets.Properties.VariableNames,'ActiveSets'));
nameAS = ActiveSets.Properties.VariableNames(contains(ActiveSets.Properties.VariableNames,'ActiveSets'));
nvarY = DSTEM_model.stem_data.stem_varset_p.nvar;
namevarY = DSTEM_model.stem_data.stem_varset_p.Y_name;
if contains(cv_part.cv_type,'LeaveOneStatOut')
    varY_pos = horzcat(cv_part.idx_st_poll{:});
else
    varY_pos = ones(size(cv_part.ARPA_stats_reg,1),nvarY);
end
MSE = cell(nvarY,cv_part.K);
RMSE = cell(nvarY,cv_part.K);
MAE = cell(nvarY,cv_part.K);
R2 = cell(nvarY,cv_part.K);
fold_failed = zeros(1,cv_part.K);
parfor run = 1:cv_part.K
    sprintf(['%%%%%%%%%%%%%%%%%%%%%%%%%%%% ',...
    'Fold %d of %d',...
    ' %%%%%%%%%%%%%%%%%%%%%%%%%%%%'],run,cv_part.K)
    %%%%% Extracting CV objects
    Y_obs = struct();
    Y_obs = cv_part.Y;
    test_coords = cv_part.test_coords{1,run};
    %%%%% STEM settings
    %%% DSTEM_varset
    % contains the observed data of all the variables and the loading coefficients;
    Covs_sel_pos = cell(1,nvarY);
    X_beta_temp = cell(1,nvarY);
    X_beta_name_temp = cell(1,nvarY);
    for p = 1:nvarY
        Covs_sel_pos{1,p} = find(ActiveSets{ActiveSets.Response_idx == p,ASi_cols(ASi_idx)});
        X_beta_temp{1,p} = DSTEM_model.stem_data.stem_varset_p.X_beta{1,p}(:,Covs_sel_pos{1,p},:);
        X_beta_name_temp{1,p} = DSTEM_model.stem_data.stem_varset_p.X_beta_name{1,p}(:,Covs_sel_pos{1,p},:);
    end
    % Y_cv - X_beta_temp*beta
    obj_stem_varset_p = stem_varset(cv_part.Y_cv{1,run},...
        DSTEM_model.stem_data.stem_varset_p.Y_name,...
        [], [], X_beta_temp, ...
        X_beta_name_temp, ...
        DSTEM_model.stem_data.stem_varset_p.X_z,...
        DSTEM_model.stem_data.stem_varset_p.X_z_name);
    %%% datestamp
    obj_temp_stem_datestamp = stem_datestamp(DSTEM_model.stem_data.stem_datestamp.date_start,...
        DSTEM_model.stem_data.stem_datestamp.date_end,DSTEM_model.T);
    %%%%% STEM_data object creation
    shape = [];
    obj_stem_modeltype = stem_modeltype(DSTEM_model.stem_data.stem_modeltype.model_name);
    obj_temp_stem_data = stem_data(obj_stem_varset_p, ...
        DSTEM_model.stem_data.stem_gridlist_p, [], [],...
        obj_temp_stem_datestamp, ...
        [], obj_stem_modeltype, shape);
    %%%%% DSTEM_par object creation
    %%% Parameter constraints
    obj_stem_par_constraints = stem_par_constraints();
    obj_stem_par_constraints.time_diagonal = 0;
    obj_temp_stem_par = stem_par(obj_temp_stem_data, 'exponential',...
        obj_stem_par_constraints);
    %%%%% DSTEM_model object creation
    obj_temp_stem_model = stem_model(obj_temp_stem_data, obj_temp_stem_par);
    %%%%% Data transform
    if DSTEM_model.stem_data.log_transformed == 1
        obj_temp_stem_model.stem_data.log_transform;
    end
    if DSTEM_model.stem_data.standardized == 1
        obj_temp_stem_model.stem_data.standardize;
    end
    %%%%% Model estimation
    %%% Starting values
    obj_temp_stem_par.beta = DSTEM_model.stem_EM_result.stem_par.beta(find(ActiveSets{:,ASi_cols(ASi_idx)}));
    obj_temp_stem_par.theta_z = DSTEM_model.stem_EM_result.stem_par.theta_z;
    obj_temp_stem_par.v_z = DSTEM_model.stem_EM_result.stem_par.v_z;
    obj_temp_stem_par.sigma_eta = DSTEM_model.stem_EM_result.stem_par.sigma_eta;
    obj_temp_stem_par.G = DSTEM_model.stem_EM_result.stem_par.G;
    obj_temp_stem_par.sigma_eps = DSTEM_model.stem_EM_result.stem_par.sigma_eps;
    obj_temp_stem_model.set_initial_values(obj_temp_stem_par);
    %%% Parameters estimation
    obj_temp_stem_EM_options = stem_EM_options();
    obj_temp_stem_EM_options.exit_tol_loglike = DSTEM_model.stem_EM_result.exit_tol_loglike;
    obj_temp_stem_EM_options.exit_tol_par = DSTEM_model.stem_EM_result.exit_tol_par;
    obj_temp_stem_EM_options.max_iterations = DSTEM_model.stem_EM_result.iterations;
    % obj_temp_stem_EM_options.max_iterations = 1;
    y_oos = struct();
    y_oos.MSE = cell(nvarY,1);
    y_oos.RMSE = cell(nvarY,1);
    y_oos.MAE = cell(nvarY,1);
    y_oos.R2 = cell(nvarY,1);
    try
        obj_temp_stem_model.EM_estimate(obj_temp_stem_EM_options);
        obj_temp_stem_model.set_logL;
        %%%%% Store CV performance indices for each run for a given Active Set
        for p = find(varY_pos(run,:))
            for i = 1:size(cv_part.test_coords{1,run}{1,p},1)
                y_oos.(namevarY{p})(i,1) = cv_part.test_coords{1,run}{1,p}(i,1);
                y_oos.(namevarY{p})(i,2) = cv_part.test_coords{1,run}{1,p}(i,2);
                y_oos.(namevarY{p})(i,3) = Y_obs{1,p}(cv_part.test_coords{1,run}{1,p}(i,1),...
                    cv_part.test_coords{1,run}{1,p}(i,2));
                y_oos.(namevarY{p})(i,4) = obj_temp_stem_model.stem_EM_result.y_hat_back{1,p}(cv_part.test_coords{1,run}{1,p}(i,1),...
                    cv_part.test_coords{1,run}{1,p}(i,2));
            end
            y_oos.(namevarY{p}) = array2table(y_oos.(namevarY{p}));
            y_oos.(namevarY{p}).Properties.VariableNames = {'row','col','y_obs','y_hat'};
            y_oos.MSE{p,1} = nanmean((y_oos.(namevarY{p}).y_hat - y_oos.(namevarY{p}).y_obs).^2);
            y_oos.RMSE{p,1} = sqrt(y_oos.MSE{p,1});
            y_oos.MAE{p,1} = nanmean(abs(y_oos.(namevarY{p}).y_hat - y_oos.(namevarY{p}).y_obs));
            y_oos.R2{p,1} = 1 - y_oos.MSE{p,1} ./ nanvar(y_oos.(namevarY{p}).y_obs);
        end
        sprintf('Fold %d complete', run)
        AIC(run) = obj_temp_stem_model.stem_EM_result.AIC;
        LogL(run) = obj_temp_stem_model.stem_EM_result.logL;
        MSE(:,run) = y_oos.MSE;
        RMSE(:,run) = y_oos.RMSE;
        MAE(:,run) = y_oos.MAE;
        R2(:,run) = y_oos.R2;
    catch
        sprintf('Fold %d failed', run)
        fold_failed(run) = 1;
        AIC(run) = NaN;
        LogL(run) = NaN;
        MSE(:,run) = y_oos.MSE;
        RMSE(:,run) = y_oos.RMSE;
        MAE(:,run) = y_oos.MAE;
        R2(:,run) = y_oos.R2;
        continue
    end
end

CV_perf_metrics.AIC = AIC;
CV_perf_metrics.LogL = LogL;
CV_perf_metrics.MSE = MSE;
CV_perf_metrics.RMSE = RMSE;
CV_perf_metrics.MAE = MAE;
CV_perf_metrics.R2 = R2;
CV_perf_metrics.Fold_failed = fold_failed;
CV_perf_metrics.NameVarY = namevarY;
CV_perf_metrics.NameAS = nameAS;

return;
end