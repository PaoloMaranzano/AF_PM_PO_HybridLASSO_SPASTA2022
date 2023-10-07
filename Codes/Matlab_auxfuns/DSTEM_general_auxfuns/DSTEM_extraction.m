%% %%%%%%%%%% DSTEM_extraction
%%% Extracts and tabulates the main information/results for an estimated DSTEM model.

function [MLE,Perform_metrics,Setup] = DSTEM_extraction(DSTEM_model,Y_std,X_std)

% DSTEM_model = obj_stem_model;
% DSTEM_model = obj_stem_model;
% DSTEM_model = M34_estim;
%%% Model setup information
Setup.n_beta = DSTEM_model.stem_EM_result.stem_par.n_beta;
if ismember(DSTEM_model.stem_par.stem_modeltype.model_name,"f-HDGM")
    Setup.n_dep_var = 1;
    Setup.Y_names = {DSTEM_model.stem_data.stem_varset_p.Y_name{1,1}};
    Setup.X_names = DSTEM_model.stem_data.stem_varset_p.X_beta_name{1,1};
    Setup.n_sites = DSTEM_model.stem_data.stem_varset_p.dim(1,1);
    Setup.Y_std = DSTEM_model.stem_data.stem_varset_p.Y_stds;
    Setup.X_std = DSTEM_model.stem_data.stem_varset_p.X_beta_stds;
    Setup.T = DSTEM_model.stem_data.stem_varset_p.T;
    Setup.H = DSTEM_model.stem_data.stem_varset_p.nvar;
    %%% Variable names for each response variable
    var_names = {DSTEM_model.stem_data.X_beta_name{1,1}};
    %%% 
    Setup.basis_per_var = length(var_names{1,1}) / length(Setup.X_names);
else
    Setup.n_dep_var = DSTEM_model.stem_data.stem_varset_p.nvar;
    Setup.Y_names = DSTEM_model.stem_data.stem_varset_p.Y_name;
    Setup.X_names = DSTEM_model.stem_data.stem_varset_p.X_beta_name;
    Setup.n_sites = DSTEM_model.stem_data.stem_varset_p.dim;
    Setup.Y_std = DSTEM_model.stem_data.stem_varset_p.Y_stds;
    Setup.X_std = DSTEM_model.stem_data.stem_varset_p.X_beta_stds;
    Setup.T = DSTEM_model.stem_data.stem_varset_p.T;
    %%% Variable names for each response variable
    var_names = DSTEM_model.stem_data.stem_varset_p.X_beta_name;
end
Variable = horzcat(var_names{1,:})';
for j = 1:Setup.n_dep_var
    beta_y{1,j} = repelem(j,length(var_names{1,j}));
end
beta_y = horzcat(beta_y{1,:})';


%%% Maximum likelihood estimates
% Covariates names
MLE.Reg_pars.Reg_names = var_names;
Reg_names_tab = table(Variable);
Reg_names_tab.Variable = categorical(Reg_names_tab.Variable);
Reg_names_tab.Dep_var_idx = beta_y;
MLE.Reg_pars.Reg_names_tab = Reg_names_tab;
cr_names = strcat("y",num2str(MLE.Reg_pars.Reg_names_tab.Dep_var_idx),...
    "_",string(MLE.Reg_pars.Reg_names_tab.Variable));
% VarCov for the regression parameters
MLE.VarCov = DSTEM_model.stem_EM_result.stem_par.varcov;
if ~isempty(MLE.VarCov)
    MLE.Reg_pars.VarCov_reg = MLE.VarCov(1:Setup.n_beta,1:Setup.n_beta);
    MLE.Reg_pars.VarCov_reg = array2table(MLE.Reg_pars.VarCov_reg,...
        'RowNames',cr_names,'VariableNames',cr_names);
end
% Estimates of the regression parameters
MLE.Reg_pars.Beta_vec = DSTEM_model.stem_EM_result.stem_par.beta;
for j = 1:Setup.n_dep_var
    %%% Coefficients (originals)
    Var = var_names{1,j}';
    b(:,j) = MLE.Reg_pars.Beta_vec(beta_y==j);
    if ~isempty(MLE.VarCov)
        SEb(:,j) = sqrt(diag(table2array(MLE.Reg_pars.VarCov_reg(find(beta_y == j),find(beta_y == j)))));
    end
    betas.(Setup.Y_names{j}) = table(Var);
    betas.(Setup.Y_names{j}).Var = categorical(betas.(Setup.Y_names{j}).Var);
    betas.(Setup.Y_names{j}).Coef = b(:,j);
    if ~isempty(MLE.VarCov)
        betas.(Setup.Y_names{j}).SE = SEb(:,j);
        betas.(Setup.Y_names{j}).abs_t_stat = abs(betas.(Setup.Y_names{j}).Coef ./ betas.(Setup.Y_names{j}).SE);
        %%% Statistical significance
        betas.(Setup.Y_names{j}).signif = repmat(".",length(var_names{1,j}),1);
        betas.(Setup.Y_names{j}).signif(betas.(Setup.Y_names{j}).abs_t_stat > norminv(0.95)) = "*";
        betas.(Setup.Y_names{j}).signif(betas.(Setup.Y_names{j}).abs_t_stat > norminv(0.975)) = "**";
        betas.(Setup.Y_names{j}).signif(betas.(Setup.Y_names{j}).abs_t_stat > norminv(0.99)) = "***";
        betas.(Setup.Y_names{j}).signif = categorical(betas.(Setup.Y_names{j}).signif);
    end
    
    if DSTEM_model.stem_data.standardized == 1
        %%% Coefficients (de-standardized)
        if ~isempty(Y_std)
            b_destd(:,j) = MLE.Reg_pars.Beta_vec(beta_y==j) * Y_std{j} ./ X_std{j}';
            SEb_destd(:,j) = sqrt(diag(table2array(MLE.Reg_pars.VarCov_reg(find(beta_y == j),find(beta_y == j))))) * Y_std{j} ./ X_std{j}';
        else
            if ismember(DSTEM_model.stem_par.stem_modeltype.model_name,"f-HDGM")
                A = repelem(Setup.Y_std{j} ./ Setup.X_std{j},Setup.basis_per_var);
                b_destd(:,j) = MLE.Reg_pars.Beta_vec(beta_y==j) .* A';
            else
                A = Setup.Y_std{j} ./ Setup.X_std{j};
                b_destd(:,j) = MLE.Reg_pars.Beta_vec(beta_y==j) * Setup.Y_std{j} ./ Setup.X_std{j}';
            end
            vc = MLE.Reg_pars.VarCov_reg(find(beta_y == j),find(beta_y == j));
            vc_destd = A' .* table2array(vc) .* A;
            SEb_destd(:,j) = sqrt(diag(vc_destd));
            % SEb_destd(:,j) = sqrt(diag(table2array(MLE.Reg_pars.VarCov_reg(find(beta_y == j),find(beta_y == j))))) * Setup.Y_std{j} ./ Setup.X_std{j}';
        end
        betas_destd.(Setup.Y_names{j}) = table(Var);
        betas_destd.(Setup.Y_names{j}).Var = categorical(betas_destd.(Setup.Y_names{j}).Var);
        betas_destd.(Setup.Y_names{j}).Coef = b_destd(:,j);
        betas_destd.(Setup.Y_names{j}).SE = SEb_destd(:,j);
        betas_destd.(Setup.Y_names{j}).abs_t_stat = abs(betas_destd.(Setup.Y_names{j}).Coef ./ betas_destd.(Setup.Y_names{j}).SE);
        %%% Statistical significance
        betas_destd.(Setup.Y_names{j}).signif = repmat(".",length(var_names{1,j}),1);
        betas_destd.(Setup.Y_names{j}).signif(betas_destd.(Setup.Y_names{j}).abs_t_stat > norminv(0.95)) = "*";
        betas_destd.(Setup.Y_names{j}).signif(betas_destd.(Setup.Y_names{j}).abs_t_stat > norminv(0.975)) = "**";
        betas_destd.(Setup.Y_names{j}).signif(betas_destd.(Setup.Y_names{j}).abs_t_stat > norminv(0.99)) = "***";
        betas_destd.(Setup.Y_names{j}).signif = categorical(betas_destd.(Setup.Y_names{j}).signif);
        %%% Destandardized varcov
        vc_destd = array2table(vc_destd);
        vc_destd.Properties.VariableNames = Var;
        vc_destd.Properties.RowNames = Var;
        MLE.Reg_pars.VarCov_reg_destd.(Setup.Y_names{j}) = vc_destd;
    end
    
    clearvars b SEb Var b_destd SEb_destd
end
MLE.Reg_pars.Beta_tab = betas;
if DSTEM_model.stem_data.standardized == 1
    MLE.Reg_pars.Beta_destd_tab = betas_destd;
end
%%% Spatio-temporal parameters
MLE.SpatTime_pars.G = DSTEM_model.stem_EM_result.stem_par.G;
MLE.SpatTime_pars.Theta_Z = DSTEM_model.stem_EM_result.stem_par.theta_z;
MLE.SpatTime_pars.V_z = DSTEM_model.stem_EM_result.stem_par.v_z;
MLE.SpatTime_pars.Sigma_eta = DSTEM_model.stem_EM_result.stem_par.sigma_eta;
MLE.SpatTime_pars.Sigma_eps = DSTEM_model.stem_EM_result.stem_par.sigma_eps;
%%% In-sample performance metrics
Perform_metrics.AIC = DSTEM_model.stem_EM_result.AIC;
Perform_metrics.LogL = DSTEM_model.stem_EM_result.logL;
Perform_metrics.R2_InSample = DSTEM_model.stem_EM_result.R2;

end

