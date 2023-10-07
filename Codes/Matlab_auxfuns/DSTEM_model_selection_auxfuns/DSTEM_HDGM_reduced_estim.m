function [DSTEM_model_estim] = DSTEM_HDGM_reduced_estim(DSTEM_model,DSTEM_setup,ActiveSets,ASi_idx)

if 0
    DSTEM_model = obj_stem_model_full;
    ASi_idx = 1;
    DSTEM_setup = Setup;
end

%%% Model settings
if ~isfield(DSTEM_setup,'log_transform')
    DSTEM_setup.log_transform = DSTEM_model.stem_data.log_transformed;
end
if ~isfield(DSTEM_setup,'standardize')
    DSTEM_setup.standardize = DSTEM_model.stem_data.standardized;
end
if ~isfield(DSTEM_setup,'exit_tol_par')
    DSTEM_setup.exit_tol_par = DSTEM_model.stem_EM_result.exit_tol_par/2;
end
if ~isfield(DSTEM_setup,'exit_tol_loglike')
    DSTEM_setup.exit_tol_loglike = DSTEM_model.stem_EM_result.exit_tol_loglike/2;
end
if ~isfield(DSTEM_setup,'max_iterations')
    DSTEM_setup.max_iterations = DSTEM_model.stem_EM_result.max_iterations+100;
end
if ~isfield(DSTEM_setup,'LogLik_compute')
    DSTEM_setup.LogLik_compute = 0;
end
if ~isfield(DSTEM_setup,'VarCov_compute')
    DSTEM_setup.VarCov_compute = 0;
end


ASi_cols = find(contains(ActiveSets.Properties.VariableNames,'ActiveSets'));
nvarY = DSTEM_model.stem_data.stem_varset_p.nvar;



%%%%% STEM settings
%%% DSTEM_varset
% contains the observed data of all the variables and the loading coefficients;
Covs_sel_pos = cell(1,nvarY);
X_beta_temp = cell(1,nvarY);
X_beta_name_temp = cell(1,nvarY);

for p = 1:nvarY
    Covs_sel_pos{1,p} = find(ActiveSets{ActiveSets.Response_idx == p,ASi_cols(ASi_idx)});
    if DSTEM_setup.standardize == 1
        for j = 1:size(DSTEM_model.stem_data.stem_varset_p.X_beta{1,p},2)
            X_beta{1,p}(:,j,:) = DSTEM_model.stem_data.stem_varset_p.X_beta{1,p}(:,j,:) * ...
                DSTEM_model.stem_data.stem_varset_p.X_beta_stds{1,p}(j) + ...
                DSTEM_model.stem_data.stem_varset_p.X_beta_means{1,p}(j);
        end
        X_beta_temp{1,p} = X_beta{1,p}(:,Covs_sel_pos{1,p},:);
    else
        X_beta_temp{1,p} = DSTEM_model.stem_data.stem_varset_p.X_beta{1,p}(:,Covs_sel_pos{1,p},:);
    end
    X_beta_name_temp{1,p} = DSTEM_model.stem_data.stem_varset_p.X_beta_name{1,p}(:,Covs_sel_pos{1,p},:);
end

if DSTEM_setup.standardize == 1
    for p = 1:nvarY
        Y_temp{1,p} = DSTEM_model.stem_data.stem_varset_p.Y{1,p} * ...
            DSTEM_model.stem_data.stem_varset_p.Y_stds{1,p} + ...
            DSTEM_model.stem_data.stem_varset_p.Y_means{1,p};
    end
else
    Y_temp = DSTEM_model.stem_data.stem_varset_p.Y;
end

obj_stem_varset_p = stem_varset(Y_temp,...
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
if DSTEM_setup.log_transform == 1
    obj_temp_stem_model.stem_data.log_transform;
end
if DSTEM_setup.standardize == 1
    obj_temp_stem_model.stem_data.standardize;
end
%%%%% Model estimation
%%% Starting values
% obj_temp_stem_par.beta = DSTEM_model.stem_EM_result.stem_par.beta(find(ActiveSets{:,ASi_cols(ASi_idx)}));
% obj_temp_stem_par.theta_z = DSTEM_model.stem_EM_result.stem_par.theta_z;
% obj_temp_stem_par.v_z = DSTEM_model.stem_EM_result.stem_par.v_z;
% obj_temp_stem_par.sigma_eta = DSTEM_model.stem_EM_result.stem_par.sigma_eta;
% obj_temp_stem_par.G = DSTEM_model.stem_EM_result.stem_par.G;
% obj_temp_stem_par.sigma_eps = DSTEM_model.stem_EM_result.stem_par.sigma_eps;

beta_init = DSTEM_model.get_beta0();
obj_temp_stem_par.beta = beta_init(find(ActiveSets{:,ASi_cols(ASi_idx)}));
obj_temp_stem_par.theta_z = DSTEM_model.stem_par_initial.theta_z;
obj_temp_stem_par.v_z = DSTEM_model.stem_par_initial.v_z;
obj_temp_stem_par.sigma_eta = DSTEM_model.stem_par_initial.sigma_eta;
obj_temp_stem_par.G = DSTEM_model.stem_par_initial.G;
obj_temp_stem_par.sigma_eps = DSTEM_model.stem_par_initial.sigma_eps;

obj_temp_stem_model.set_initial_values(obj_temp_stem_par);

%%% Parameters estimation
obj_temp_stem_EM_options = stem_EM_options();
obj_temp_stem_EM_options.exit_tol_loglike = DSTEM_setup.exit_tol_loglike;
obj_temp_stem_EM_options.exit_tol_par = DSTEM_setup.exit_tol_par;
obj_temp_stem_EM_options.max_iterations = DSTEM_setup.max_iterations;
obj_temp_stem_model.EM_estimate(obj_temp_stem_EM_options);
%%% Log-likelihood computation
if DSTEM_setup.LogLik_compute == 1
    obj_temp_stem_model.set_logL;
end
%%% VarCov matrix computation
if DSTEM_setup.VarCov_compute == 1
    obj_temp_stem_model.set_varcov;
end


DSTEM_model_estim = obj_temp_stem_model;



return;
end