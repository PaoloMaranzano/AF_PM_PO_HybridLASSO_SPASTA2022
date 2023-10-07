function [ActiveSets,nAS,beta_LASSO,beta_LASSO_tab,FitInfo_LASSO] = DSTEM_SepLASSO(...
    DSTEM_model_full,Lambda_seq,True_model)

debug = 0;
if debug == 1
   DSTEM_model_full = obj_stem_model_full; 
end
sprintf(['Separated LASSO Estimation begin at ' datestr(datetime('now'))])
%%% Response variables
nY =length(DSTEM_model_full.stem_data.stem_varset_p.Y);
%%% Covariates names
var_names = DSTEM_model_full.stem_data.stem_varset_p.X_beta_name;
%%% Covariates indices
if debug == 1
    idx_betas = ismember(horzcat(var_names{1,:}),'Constant');
else
    idx_betas = ismember(horzcat(var_names{1,:}),'cons');
end
%%% Covariates positions
pos_betas = find(~idx_betas);
%%% Covariates names
X_names = horzcat(var_names{1,:});
%%% Intercepts positions
pos_const = find(idx_betas);
%%% ML estimates of the regression parameters
beta_MLE = DSTEM_model_full.stem_EM_result.stem_par.beta(pos_betas);



%% %%%%% 1-dim matrices building
for p = 1:nY
    %%% Dependent variable
    % Original observations (2-dimensions)
    y_2d{p} = DSTEM_model_full.stem_data.stem_varset_p.Y{p};
    ncols = size(y_2d{p},2);
    nrows = size(y_2d{p},1);
    % Reshape to 1-dim (stacked obs)
    y_1d{p} = reshape(y_2d{p}.',1,[])';
    %%% Covariates
    % Original observations (3-dimensions)
    X_3d{p} = DSTEM_model_full.stem_data.stem_varset_p.X_beta{p}(:,2:end,:);
    n_covs = size(X_3d{p},2);
    clearvars Xv_1d Xv_2d
    for v = 1:n_covs
        Xv_2d = reshape(X_3d{p}(:,v,:),nrows,ncols,[]);
        Xv_1d(:,v) = reshape(Xv_2d.',1,[])';
    end
    X_1d{p} = Xv_1d;
    % g{p} = repelem(p,n_covs);
end
% Stack Y
y = vertcat(y_1d{:});
% Block-diagonalization of the covariates
X = blkdiag(X_1d{:});
% Stack indices
% index = horzcat(g{:});
% Drop missing values
nan_pos = isnan(y);
y = y(~nan_pos);
X = X(~nan_pos,:);



%% %%%%% Separated LASSO
LASSO_begin_time = datetime('now')
parfor lam = 1:length(Lambda_seq)
    sprintf('Lambda %d of %d', lam, length(Lambda_seq))
    beta_LASSO_p = cell(1,nY);
    FitInfo_p = cell(1,nY);
    for p = 1:nY
        % 'CV','resubstitution'
        % [B,FitInfo] = lasso(X_1d{p},y_1d{p},'CV','resubstitution',...
        %     'Standardize',true,'MaxIter',1e5,...
        %     'Options',statset('UseParallel',true));
        'CV','resubstitution'
        [B,FitInfo] = lasso(X_1d{p},y_1d{p},'CV',10,...
            'Lambda',Lambda_seq(lam),...
            'Standardize',true,'MaxIter',1e5,...
            'Options',statset('UseParallel',false));
        beta_LASSO_p{p} = B;
        FitInfo_p{p} = FitInfo;
    end
   beta_LASSO_lam = vertcat(beta_LASSO_p{:});
   beta_LASSO(:,lam) = beta_LASSO_lam;
   FitInfo_LASSO{lam} = FitInfo_p;
end
beta_LASSO_tab = table(X_names(pos_betas)');
beta_LASSO_tab = [beta_LASSO_tab , array2table(beta_LASSO)];
for lam = 1 : length(Lambda_seq)
    lam_num{lam} = ['Lambda_' num2str(Lambda_seq(lam),3)];        % Character Array
end
beta_LASSO_tab.Properties.VariableNames = {'Variable',lam_num{:}}
LASSO_end_time = datetime('now')



%% %%%%% ActiveSets definition
sprintf('Active Sets definition - Begin')
%%% Find coefficients above the treshold (non null)
ActiveCovs = beta_LASSO ~= 0;
%%% Define the Active Sets
ActiveSets_temp = unique(ActiveCovs', 'rows')';
nAS = size(ActiveSets_temp,2);
for p = 1:nY
    if debug == 1
        beta_y{1,p} = repelem(p,sum(~ismember(var_names{1,p},'Constant')));
    else
        beta_y{1,p} = repelem(p,sum(~ismember(var_names{1,p},'cons')));
    end    
end
beta_y = horzcat(beta_y{1,:})';
%%% Adding intercept to the ActiveSets and beta_y
for ASi = 1:nAS
    for p = 1:nY
        ActiveCovs_y{1,p} = [1,ActiveSets_temp(find(beta_y == p),ASi)'];
        beta_y_full{1,p} = repelem(p,length(var_names{1,p}));
    end
    ActiveSets(:,ASi) = horzcat(ActiveCovs_y{:})';
end
%%% Convert ActiveSets to table
Variable = horzcat(var_names{1,:})';
Response_idx = horzcat(beta_y_full{:})';
if isempty(True_model) == 0
    for ASi = 1:nAS
        true_model_check(1,ASi) = all(ActiveSets(:,ASi) == True_model);
    end
    if find(true_model_check) ~= 0
        AS_true = ['True model equals to ActiveSets' num2str(find(true_model_check))];
    else
        AS_true = 'True model not in the ActiveSets';
    end
    ActiveSets = [table(Variable), ...
        table(Response_idx), table(True_model), array2table(ActiveSets)];
else
    ActiveSets = [table(Variable), ...
        table(Response_idx), array2table(ActiveSets)];
end
sprintf('Active Sets definition - End')
sprintf(['Separated LASSO Estimation ended at ' datestr(datetime('now'))])

    return;
end