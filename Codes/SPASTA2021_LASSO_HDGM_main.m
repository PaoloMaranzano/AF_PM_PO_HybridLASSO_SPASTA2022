if 1
    clear
    clc
    addpath(genpath('../../../../SPASTA2021/'));
    % load('Output/Output_09052021/CV_part.mat')
    % load('Output/Output_10052021/Full_model.mat')
    % load('Output/Output_09052021/LASSO.mat')
    % load('Output/Output_09052021/LOSO_CV.mat')
    % load('Output/Output_09052021/randK_CV.mat')
    % load('Output/Output_09052021/stratK_CV.mat')
    clearvars data_preparation standardization log_transform save_all save_parts full_estimation
    clearvars full_varcov lasso_flag parfor_ext_flag parfor_int_1loop_flag parfor_int_3loop_flag
    clearvars randK_cv_flag stratK_cv_flag loso_cv_flag
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%% HDGM with LASSO model seleciton %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Import and reshape AQ data for Lombardy 
if ~exist('Ground')
    addpath(genpath('../../../../SPASTA2021/'));
    run('Reshape_to_HDGM.m');
end
% clearvars -except crossval_step data Ground log_transform poll standardization ActiveSets nAS beta_LASSO beta_LASSO_tab FitInfo_LASSO

% Settings
data_preparation = 1;
standardization = 1;
log_transform = 0;
save_parts = 1;
save_all = 1;
save_timing = 1;
full_estimation = 1;
full_varcov = 1;
lasso_flag = 1;
grid_dense = 0;

parfor_int_1loop_flag = 0;
parfor_int_3loop_flag = 0;
randK_cv_parfor_int_flag = 0;
stratK_cv_parfor_int_flag = 0;
loso_cv_parfor_int_flag = 0;

parfor_ext_flag = 1;
randK_cv_parfor_ext_flag = 1;
stratK_cv_parfor_ext_flag = 0;
loso_cv_parfor_ext_flag = 0;

time_begin = datetime('now')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%      Step 0. Data definition     %%
%%    ECMWF + ARPA Lombardia data   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if data_preparation == 1
    sprintf('Step 0: data definition')
    time_data_definition_begin = datetime('now')
    %% Data building
    poll = Ground.poll;
    n = cell(length(poll),1);
    obj_stem_gridlist_p = stem_gridlist();
    for p = 1:length(poll)
        %%% Dependent variables
        ground.Y{p} = Ground.([poll{p}]);
        ground.Y_name{p} = [poll{p} '_ground'];
        n{p,1} = size(ground.Y{p}, 1);
        %%% Loading covariates for the selected pollutants at each monitoring stations
        ground.X_beta{p} = Ground.(['X_' poll{p}]);
        ground.X_beta_name{p} = Ground.vars_names;
        ground.X_z{p} = ones(n{p}, 1);
        ground.X_z_name{p} = {['constant_' poll{p}]};
        %%% Coordinates grid
        ground.coordinates{p} = [Ground.(['Coords_' poll{p}]).Latitude, ...
            Ground.(['Coords_' poll{p}]).Longitude];
        obj_stem_grid = cell(length(poll),1);
        obj_stem_grid{p} = stem_grid(ground.coordinates{p}, 'deg', 'sparse', 'point');
        %%% stem_gridlist
        obj_stem_gridlist_p.add(obj_stem_grid{p});
    end
    T = size(ground.Y{1}, 2);
    %% DSTEM settings
    %%% DSTEM_varset
    obj_stem_varset_p = stem_varset(ground.Y, ground.Y_name, [], [], ...
        ground.X_beta, ground.X_beta_name, ...
        ground.X_z, ground.X_z_name);
    %%% DSTEM_data
    obj_stem_datestamp = stem_datestamp('01-01-2017 00:00','31-10-2020 00:00',T);
    shape = [];
    obj_stem_modeltype = stem_modeltype('HDGM');
    obj_stem_data = stem_data(obj_stem_varset_p, obj_stem_gridlist_p, ...
        [], [], obj_stem_datestamp, [], ...
        obj_stem_modeltype, shape);
    %%% DSTEM_par
    obj_stem_par_constr_red=stem_par_constraints();
    obj_stem_par_constr_red.time_diagonal=0;
    obj_stem_par = stem_par(obj_stem_data, 'exponential',obj_stem_par_constr_red);
    obj_stem_model = stem_model(obj_stem_data, obj_stem_par);
    %%% Data transform
    if log_transform == 1
        obj_stem_model.stem_data.log_transform;
    end
    if standardization == 1
        obj_stem_model.stem_data.standardize;
    end
    time_data_definition_end = datetime('now')
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step 1. Full model estimation %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if full_estimation == 1
    sprintf('Step 1: full model estimation')
    obj_stem_model_full = obj_stem_model;
    %%% Starting values
    obj_stem_par.beta = obj_stem_model_full.get_beta0();
    obj_stem_par.theta_z = km2deg(50);         % Kilometers
    obj_stem_par.v_z = [
        1 0.6 0.6;
        0.6 1 0.9
        0.6 0.9 1];             % Cross-correlation between multiple variables
    obj_stem_par.sigma_eta = diag([0.2 0.2 0.2]);
    obj_stem_par.G = diag([0.8 0.8 0.8]);
    obj_stem_par.sigma_eps = diag([0.3 0.3 0.3]);
    obj_stem_model_full.set_initial_values(obj_stem_par);
    %%% Parameters estimation
    obj_stem_EM_options = stem_EM_options();
    obj_stem_EM_options.exit_tol_loglike = 0.0001;
    obj_stem_EM_options.exit_tol_par = 0.0001;
    obj_stem_EM_options.max_iterations = 200;
    time_full_estim_begin = datetime('now')
    obj_stem_model_full.EM_estimate(obj_stem_EM_options);
    time_full_estim_end = datetime('now')
    if full_varcov == 1
        %%% Exact VarCov matrix of estimated pars
        time_full_varcov_begin = datetime('now')
        obj_stem_model_full.set_varcov;
        time_full_varcov_end = datetime('now')
        %%% Log-likelihood of the data
        time_full_loglik_begin = datetime('now')
        obj_stem_model_full.set_logL;
        time_full_loglik_end = datetime('now')
        %%% Extract main objects from full model
        [Full_model_MLE,Full_model_pef_met,...
            Full_model_setup] = DSTEM_extraction(obj_stem_model_full);
    end
    if save_parts == 1
        save('Full_model.mat')
        % save('Full_model_NonStd.mat')
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step 2. LASSO and Active Sets definition %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if lasso_flag == 1
    if grid_dense == 0
        Lambda_seq = [0 , [0.0001:0.0002:1], [1.01:0.01:4]];
        % Lambda_seq = [0 , [0.0001:0.0001:0.5], [0.51:0.01:4]];
    end
    if grid_dense == 1
        Lambda_seq = linspace(0,0.0001*4,50);
    end
    sprintf('Step 3a: Separated LASSO')
    time_LASSO_begin = datetime('now')
    [ActiveSets,nAS,beta_LASSO,beta_LASSO_tab,FitInfo_LASSO] = DSTEM_SepLASSO(...
        obj_stem_model_full,Lambda_seq,[]);
    time_LASSO_end = datetime('now')
    if save_parts == 1
        save('LASSO.mat','ActiveSets','nAS','beta_LASSO','beta_LASSO_tab',...
            'FitInfo_LASSO','Lambda_seq','time_LASSO_begin','time_LASSO_end')
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step 3. Cross-validation partitioning %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sprintf('Step 2: cross-validation partitioning')
time_CVpart_begin = datetime('now')
%%% Number of folds
K = 10;
%%% Partitioning data
% Stratified K-fold CV algorithm
cv_part_stratK = CVpart_strat_Kfold(Ground,'Tipology_rec',K);
% Random K-fold CV algorithm
cv_part_randK = CVpart_random_Kfold(Ground,K);
% Leave-one-station-out CV algorithm
cv_part_LOSO = CVpart_LeaveOneStatOut(Ground);
time_CVpart_end = datetime('now')
if save_parts == 1
    save('CV_part.mat','cv_part_stratK','cv_part_randK','cv_part_LOSO',...
        'time_CVpart_begin','time_CVpart_end')
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step 4. Spatio-temporal stratified K-fold CV %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Parfor esterno
if parfor_ext_flag == 1
    time_CV_begin = datetime('now')
    if stratK_cv_parfor_ext_flag == 1
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Step 4a. Spatio-temporal stratified K-fold CV %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        sprintf('Step 4a: spatio-temporal stratified K-fold CV for each active set')
        time_StratKfoldCV_begin = datetime('now')
        parfor ASi = 1:nAS
            sprintf('Estimating Active Set number %d of %d',ASi,nAS)
            %%%%% Cross-validation
            Perf_metrics = DSTEM_KfoldCV_par(cv_part_stratK,obj_stem_model_full,ActiveSets,ASi);
            %%%%% Store CV performance indices for each of the Active Sets
            AIC_ASi(:,ASi) = Perf_metrics.AIC;
            LogL_ASi(:,ASi) = Perf_metrics.LogL;
            MSE_ASi{ASi} = Perf_metrics.MSE';
            RMSE_ASi{ASi} = Perf_metrics.RMSE';
            MAE_ASi{ASi} = Perf_metrics.MAE';
            R2_ASi{ASi} = Perf_metrics.R2';
            Fold_failed(:,ASi) = Perf_metrics.Fold_failed';
        end
        stratK_CV_metrics.AIC_ASi = AIC_ASi;
        stratK_CV_metrics.LogLASi = LogL_ASi;
        for ASi = 1:nAS
            for p = 1:length(poll)
                stratK_CV_metrics.MSE_ASi{p}(:,ASi) = MSE_ASi{ASi}(:,p);
                stratK_CV_metrics.RMSE_ASi{p}(:,ASi) = RMSE_ASi{ASi}(:,p);
                stratK_CV_metrics.MAE_ASi{p}(:,ASi) = MAE_ASi{ASi}(:,p);
                stratK_CV_metrics.R2_ASi{p}(:,ASi) = R2_ASi{ASi}(:,p);
            end
        end
        stratK_CV_metrics.Fold_failed = Fold_failed;
        time_StratKfoldCV_end = datetime('now')
        clear MSE_ASi RMSE_ASi MAE_ASi R2_ASi LogL_ASi AIC_ASi Fold_failed Perf_metrics;
        if save_parts == 1
            save('stratK_CV.mat','stratK_CV_metrics')
        end
    end
    
    if randK_cv_parfor_ext_flag == 1
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Step 4b. Spatio-temporal random K-fold CV %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        sprintf('Step 4b: spatio-temporal random K-fold CV for each active set')
        time_RandKfoldCV_begin = datetime('now')
        parfor ASi = 1:nAS
            sprintf('Estimating Active Set number %d of %d',ASi,nAS)
            %%%%% Cross-validation
            Perf_metrics = DSTEM_KfoldCV_par(cv_part_randK,obj_stem_model_full,ActiveSets,ASi);
            %%%%% Store CV performance indices for each of the Active Sets
            AIC_ASi(:,ASi) = Perf_metrics.AIC;
            LogL_ASi(:,ASi) = Perf_metrics.LogL;
            MSE_ASi{ASi} = Perf_metrics.MSE';
            RMSE_ASi{ASi} = Perf_metrics.RMSE';
            MAE_ASi{ASi} = Perf_metrics.MAE';
            R2_ASi{ASi} = Perf_metrics.R2';
            Fold_failed(:,ASi) = Perf_metrics.Fold_failed';
        end
        randK_CV_metrics.AIC_ASi = AIC_ASi;
        randK_CV_metrics.LogLASi = LogL_ASi;
        for ASi = 1:nAS
            for p = 1:length(poll)
                randK_CV_metrics.MSE_ASi{p}(:,ASi) = MSE_ASi{ASi}(:,p);
                randK_CV_metrics.RMSE_ASi{p}(:,ASi) = RMSE_ASi{ASi}(:,p);
                randK_CV_metrics.MAE_ASi{p}(:,ASi) = MAE_ASi{ASi}(:,p);
                randK_CV_metrics.R2_ASi{p}(:,ASi) = R2_ASi{ASi}(:,p);
            end
        end
        randK_CV_metrics.Fold_failed = Fold_failed;
        time_RandKfoldCV_end = datetime('now')
        clear MSE_ASi RMSE_ASi MAE_ASi R2_ASi LogL_ASi AIC_ASi Fold_failed Perf_metrics;
        if save_parts == 1
            save('randK_CV.mat','randK_CV_metrics')
        end
    end
    
    if loso_cv_parfor_ext_flag == 1
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Step 4c. Spatial Leave-one-stat-out CV %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        sprintf('Step 4c: spatial leave-one-stat-out CV for each active set')
        time_LOSO_CV_begin = datetime('now')
        parfor ASi = 1:nAS
            sprintf('Estimating Active Set number %d of %d',ASi,nAS)
            %%%%% Cross-validation
            Perf_metrics = DSTEM_KfoldCV_par(cv_part_LOSO,obj_stem_model_full,ActiveSets,ASi);
            %%%%% Store CV performance indices for each of the Active Sets
            AIC_ASi(:,ASi) = Perf_metrics.AIC;
            LogL_ASi(:,ASi) = Perf_metrics.LogL;
            MSE_ASi{ASi} = Perf_metrics.MSE';
            RMSE_ASi{ASi} = Perf_metrics.RMSE';
            MAE_ASi{ASi} = Perf_metrics.MAE';
            R2_ASi{ASi} = Perf_metrics.R2';
            Fold_failed(:,ASi) = Perf_metrics.Fold_failed';
        end
        LOSO_CV_metrics.AIC_ASi = AIC_ASi;
        LOSO_CV_metrics.LogLASi = LogL_ASi;
        for ASi = 1:nAS
            for p = 1:length(poll)
                LOSO_CV_metrics.MSE_ASi{p}(:,ASi) = MSE_ASi{ASi}(:,p);
                LOSO_CV_metrics.RMSE_ASi{p}(:,ASi) = RMSE_ASi{ASi}(:,p);
                LOSO_CV_metrics.MAE_ASi{p}(:,ASi) = MAE_ASi{ASi}(:,p);
                LOSO_CV_metrics.R2_ASi{p}(:,ASi) = R2_ASi{ASi}(:,p);
            end
        end
        LOSO_CV_metrics.Fold_failed = Fold_failed;
        time_LOSO_CV_end = datetime('now')
        clear MSE_ASi RMSE_ASi MAE_ASi R2_ASi LogL_ASi AIC_ASi Fold_failed Perf_metrics;
        if save_parts == 1
            save('LOSO_CV.mat','LOSO_CV_metrics')
        end
    end
    time_CV_end = datetime('now');
end

%%%%% Parfor interno con tre loop diversi per ogni tipologia di CV
if parfor_int_3loop_flag == 1
    if stratK_cv_parfor_int_flag == 1
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Step 4a. Spatio-temporal stratified K-fold CV %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        sprintf('Step 4a: spatio-temporal stratified K-fold CV for each active set')
        time_StratKfoldCV_begin = datetime('now')
        clear MSE_lambdas RMSE_lambdas MAE_lambdas R2_lambdas
        for ASi = 1:nAS
            sprintf('Estimating Active Set number %d of %d', ASi,nAS)
            Perf_CV_metrics = DSTEM_KfoldCV_par(cv_part_stratK,obj_stem_model_full,ActiveSets,ASi);
            %%%%% Store CV performance indices for each of the Active Sets
            stratK_CV_metrics.AIC_ASi(:,ASi) = Perf_CV_metrics.AIC;
            stratK_CV_metrics.LogL_ASi(:,ASi) = Perf_CV_metrics.LogL;
            for p = 1:length(poll)
                stratK_CV_metrics.MSE_ASi{p}(:,ASi) = Perf_CV_metrics.MSE(p,:)';
                stratK_CV_metrics.RMSE_ASi{p}(:,ASi) = Perf_CV_metrics.RMSE(p,:)';
                stratK_CV_metrics.MAE_ASi{p}(:,ASi) = Perf_CV_metrics.MAE(p,:)';
                stratK_CV_metrics.R2_ASi{p}(:,ASi) = Perf_CV_metrics.R2(p,:)';
            end
            stratK_CV_metrics.Fold_failed(:,ASi) = Perf_CV_metrics.Fold_failed';
        end
        time_StratKfoldCV_end = datetime('now')
        clear MSE RMSE MAE R2 LogL AIC Fold_failed Perf_CV_metrics;
        if save_parts == 1
            save('stratK_CV.mat','stratK_CV_metrics')
        end
    end
    
    if randK_cv_parfor_int_flag == 1
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Step 4b. Spatio-temporal random K-fold CV %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        sprintf('Step 4b: spatio-temporal random K-fold CV for each active set')
        time_RandKfoldCV_begin = datetime('now')
        clear MSE_lambdas RMSE_lambdas MAE_lambdas R2_lambdas
        for ASi = 1:nAS
            sprintf('Estimating Active Set number %d of %d', ASi,nAS)
            Perf_CV_metrics = DSTEM_KfoldCV_par(cv_part_randK,obj_stem_model_full,ActiveSets,ASi);
            %%%%% Store CV performance indices for each of the Active Sets
            randK_CV_metrics.AIC_ASi(:,ASi) = Perf_CV_metrics.AIC;
            randK_CV_metrics.LogL_ASi(:,ASi) = Perf_CV_metrics.LogL;
            for p = 1:length(poll)
                randK_CV_metrics.MSE_ASi{p}(:,ASi) = Perf_CV_metrics.MSE(p,:)';
                randK_CV_metrics.RMSE_ASi{p}(:,ASi) = Perf_CV_metrics.RMSE(p,:)';
                randK_CV_metrics.MAE_ASi{p}(:,ASi) = Perf_CV_metrics.MAE(p,:)';
                randK_CV_metrics.R2_ASi{p}(:,ASi) = Perf_CV_metrics.R2(p,:)';
            end
            randK_CV_metrics.Fold_failed(:,ASi) = Perf_CV_metrics.Fold_failed';
        end
        time_RandKfoldCV_end = datetime('now')
        clear MSE RMSE MAE R2 LogL AIC Fold_failed Perf_CV_metrics;
        if save_parts == 1
            save('randK_CV.mat','randK_CV_metrics')
        end
    end
    
    if loso_cv_parfor_int_flag == 1
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Step 4c. Spatial Leave-one-stat-out CV %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        sprintf('Step 4c: spatial leave-one-stat-out CV for each active set')
        time_LOSO_CV_begin = datetime('now')
        clear MSE_lambdas RMSE_lambdas MAE_lambdas R2_lambdas
        for ASi = 1:nAS
            sprintf('Estimating Active Set number %d of %d', ASi,nAS)
            Perf_CV_metrics = DSTEM_KfoldCV_par(cv_part_LOSO,obj_stem_model_full,ActiveSets,ASi);
            %%%%% Store CV performance indices for each of the Active Sets
            LOSO_CV_metrics.AIC_ASi(:,ASi) = Perf_CV_metrics.AIC;
            LOSO_CV_metrics.LogL_ASi(:,ASi) = Perf_CV_metrics.LogL;
            for p = 1:length(poll)
                LOSO_CV_metrics.MSE_ASi{p}(:,ASi) = Perf_CV_metrics.MSE(p,:)';
                LOSO_CV_metrics.RMSE_ASi{p}(:,ASi) = Perf_CV_metrics.RMSE(p,:)';
                LOSO_CV_metrics.MAE_ASi{p}(:,ASi) = Perf_CV_metrics.MAE(p,:)';
                LOSO_CV_metrics.R2_ASi{p}(:,ASi) = Perf_CV_metrics.R2(p,:)';
            end
            LOSO_CV_metrics.Fold_failed(:,ASi) = Perf_CV_metrics.Fold_failed';
        end
        time_LOSO_CV_end = datetime('now')
        clear MSE RMSE MAE R2 LogL AIC Fold_failed Perf_CV_metrics;
        if save_parts == 1
            save('LOSO_CV.mat','LOSO_CV_metrics')
        end
    end
end

%%%%% Parfor interno con un loop solo per tuttr le tipologie di CV
if parfor_int_1loop_flag == 1
    time_CV_begin = datetime('now')
    clear MSE_lambdas RMSE_lambdas MAE_lambdas R2_lambdas
    for ASi = 1:nAS
        sprintf('Estimating Active Set number %d of %d', ASi,nAS)
        %%% Spatio-temporal stratified K-fold CV
        if 1
            sprintf('Step 4a: spatio-temporal stratified K-fold CV for each active set')
            Perf_CV_metrics = DSTEM_KfoldCV_par(cv_part_stratK,obj_stem_model_full,ActiveSets,ASi);
            stratK_CV_metrics.AIC_ASi(:,ASi) = Perf_CV_metrics.AIC;
            stratK_CV_metrics.LogL_ASi(:,ASi) = Perf_CV_metrics.LogL;
            for p = 1:length(poll)
                stratK_CV_metrics.MSE_ASi{p}(:,ASi) = Perf_CV_metrics.MSE(p,:)';
                stratK_CV_metrics.RMSE_ASi{p}(:,ASi) = Perf_CV_metrics.RMSE(p,:)';
                stratK_CV_metrics.MAE_ASi{p}(:,ASi) = Perf_CV_metrics.MAE(p,:)';
                stratK_CV_metrics.R2_ASi{p}(:,ASi) = Perf_CV_metrics.R2(p,:)';
            end
        end
        %%% Spatio-temporal random K-fold CV
        sprintf('Step 4b: spatio-temporal random K-fold CV for each active set')
        if 1
            Perf_CV_metrics = DSTEM_KfoldCV_par(cv_part_randK,obj_stem_model_full,ActiveSets,ASi);
            randK_CV_metrics.AIC_ASi(:,ASi) = Perf_CV_metrics.AIC;
            randK_CV_metrics.LogL_ASi(:,ASi) = Perf_CV_metrics.LogL;
            for p = 1:length(poll)
                randK_CV_metrics.MSE_ASi{p}(:,ASi) = Perf_CV_metrics.MSE(p,:)';
                randK_CV_metrics.RMSE_ASi{p}(:,ASi) = Perf_CV_metrics.RMSE(p,:)';
                randK_CV_metrics.MAE_ASi{p}(:,ASi) = Perf_CV_metrics.MAE(p,:)';
                randK_CV_metrics.R2_ASi{p}(:,ASi) = Perf_CV_metrics.R2(p,:)';
            end
        end
        %%% Leave-one-stat-out CV
        if 1
            sprintf('Step 4c: spatial leave-one-stat-out CV for each active set')
            Perf_CV_metrics = DSTEM_KfoldCV_par(cv_part_LOSO,obj_stem_model_full,ActiveSets,ASi);
            LOSO_CV_metrics.AIC_ASi(:,ASi) = Perf_CV_metrics.AIC;
            LOSO_CV_metrics.LogL_ASi(:,ASi) = Perf_CV_metrics.LogL;
            for p = 1:length(poll)
                LOSO_CV_metrics.MSE_ASi{p}(:,ASi) = Perf_CV_metrics.MSE(p,:)';
                LOSO_CV_metrics.RMSE_ASi{p}(:,ASi) = Perf_CV_metrics.RMSE(p,:)';
                LOSO_CV_metrics.MAE_ASi{p}(:,ASi) = Perf_CV_metrics.MAE(p,:)';
                LOSO_CV_metrics.R2_ASi{p}(:,ASi) = Perf_CV_metrics.R2(p,:)';
            end
        end
    end
    time_CV_end = datetime('now')
end

%%%% Fine 
time_end = datetime('now');

%%%% Saving timing
if save_timing == 1
    save('timing.mat','time_*')
end



if 1
    time_begin
    time_data_definition_begin
    time_data_definition_end
    time_full_varcov_begin
    time_full_varcov_end
    time_full_loglik_begin
    time_full_loglik_end
    time_LASSO_begin
    time_LASSO_end
    % time_StratKfoldCV_begin
    % time_StratKfoldCV_end
    time_RandKfoldCV_begin
    time_RandKfoldCV_end
    % time_LOSO_CV_begin
    % time_LOSO_CV_end
    time_end
    nAS
end

if save_all == 1
    save('output_PMLE_stream.mat')
end




%% optimal Active Set definition
%%%%% Find the optimal Active Set for each pollutant
if 0
    sprintf('Step 4b: optimal active set definition')
    for p = 1:length(poll)
        opt_ASi{1,p} = find(mean(cell2mat(randK_CV_metrics.MSE_ASi{1,p}),1)' == min(mean(cell2mat(randK_CV_metrics.MSE_ASi{1,p}),1)'))
    end
    plot(mean(cell2mat(randK_CV_metrics.MSE_ASi{1,p}),1)')
end



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Step 5. Compare spatio-temporal CV metrics %%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%% Find the optimal Active Set for each pollutant
% sprintf('Step 4b: optimal active set definition')
% for p = 1:length(poll)
%     opt_ASi{1,p} = find(mean(cell2mat(MSE_ASi{1,p}),1)' == min(mean(cell2mat(MSE_ASi{1,p}),1)'))
% end
% format long
% find(0.5*mean(cell2mat(MAE_ASi{1,1}),1)' + 0.5*mean(cell2mat(MAE_ASi{1,2}),1)' == min(0.5*mean(cell2mat(MAE_ASi{1,1}),1)' + 0.5*mean(cell2mat(MAE_ASi{1,2}),1)'))
% figure
% plot(mean(cell2mat(RMSE_ASi{1,1}),1)','p')
% hold on
% for p = 2:length(poll)
%     plot(mean(cell2mat(RMSE_ASi{1,p}),1)','p')
% end
% Sel_opt_ASi = 1;
% if all(cell2mat(opt_ASi))
%     Sel_opt_ASi = mean(cell2mat(opt_ASi));
% end
% ActiveSets_opt = ActiveSets(:,[1,2,3,Sel_opt_ASi+3])
% ActiveSets_opt = ActiveSets(:,[1,2,3,7+3]);


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Step 7. Final (selected) model estimation %%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sprintf('Step 7: final (selected) model estimation')
% %%% Computation begin
% date_begin = datetime('now')
% %%% STEM settings
% Covs_sel_final_pos = cell(1,length(DSTEM_str_ground.poll));
% X_beta_final = cell(1,length(DSTEM_str_ground.poll));
% X_beta_name_final = cell(1,length(DSTEM_str_ground.poll));
% for p = 1:length(DSTEM_str_ground.poll)
%     Covs_sel_final_pos{1,p} = find(ActiveSets_opt{poll_idx == p,4});
%     X_beta_final{1,p} = obj_stem_model_full.stem_data.stem_varset_p.X_beta{1,p}(:,Covs_sel_final_pos{1,p},:);
%     X_beta_name_final{1,p} = obj_stem_model_full.stem_data.stem_varset_p.X_beta_name{1,p}(:,Covs_sel_final_pos{1,p},:);
% end
% obj_final_stem_varset_p = stem_varset(cv_part.Y,ground.Y_name,...
%     [],[],X_beta_final,X_beta_name_final,ground.X_z, ground.X_z_name);
% %%% datestamp
% obj_final_stem_datestamp = stem_datestamp(ground.datestamp_begin,...
%     ground.datestamp_end,T);
% %%% HDGM STEM_data object
% shape = [];
% obj_final_stem_modeltype = stem_modeltype('HDGM');
% obj_final_stem_data = stem_data(obj_final_stem_varset_p, ...
%     ground.obj_stem_gridlist_p, [], [], obj_final_stem_datestamp, ...
%     [], obj_final_stem_modeltype, shape);
% %%% DSTEM_par object creation
% % contains the structure and the values of the model parameters;
% obj_final_stem_par_constraints=stem_par_constraints();
% %%% time_diagonal
% % used to specify if the matrices G and Sigma_eta are diagonal or not
% obj_final_stem_par_constraints.time_diagonal=0;
% obj_final_stem_par = stem_par(obj_final_stem_data, 'exponential',...
%     obj_final_stem_par_constraints);
% %%% DSTEM_model object creation
% obj_final_stem_model = stem_model(obj_final_stem_data, obj_final_stem_par);
% %%% Data transform
% if log_transform == 1
%     obj_final_stem_model.stem_data.log_transform;
% end
% if standardization == 1
%     obj_final_stem_model.stem_data.standardize;
% end
% %%% Starting values
% idx_beta_sel_opt_ASi = find(ActiveSets_opt{:,4});
% obj_final_stem_par.beta = DSTEM_obj_sim.stem_model.stem_EM_result.stem_par.beta(idx_beta_sel_opt_ASi);
% obj_final_stem_par.theta_z = DSTEM_obj_sim.stem_model.stem_EM_result.stem_par.theta_z;
% obj_final_stem_par.v_z = DSTEM_obj_sim.stem_model.stem_EM_result.stem_par.v_z;
% obj_final_stem_par.sigma_eta = DSTEM_obj_sim.stem_model.stem_EM_result.stem_par.sigma_eta;
% obj_final_stem_par.G = DSTEM_obj_sim.stem_model.stem_EM_result.stem_par.G;
% obj_final_stem_par.sigma_eps = DSTEM_obj_sim.stem_model.stem_EM_result.stem_par.sigma_eps;
% obj_final_stem_model.set_initial_values(obj_final_stem_par);
% %%% Parameters estimation
% obj_final_stem_EM_options = stem_EM_options();
% obj_final_stem_EM_options.exit_tol_loglike = DSTEM_obj_sim.stem_model.stem_EM_result.exit_tol_loglike;
% obj_final_stem_EM_options.exit_tol_par = DSTEM_obj_sim.stem_model.stem_EM_result.exit_tol_par;
% obj_final_stem_EM_options.max_iterations = DSTEM_obj_sim.stem_model.stem_EM_result.iterations+20;
% date_estimation_begin = datetime('now')
% obj_final_stem_model.EM_estimate(obj_final_stem_EM_options);
% date_estimation_end = datetime('now')


%%%% Store for each iteration
% lambda sequencefigure(1)
% plot(mean(cell2mat(LOSO_CV_metrics.MAE_ASi{1,1}),1))
% figure(2)
% plot(mean(cell2mat(LOSO_CV_metrics.RMSE_ASi{1,1}),1))
% clearvars -except MSE_ASi RMSE_ASi MAE_ASi R2_ASi LogL_ASi AIC_ASi ActiveSets opt_ASi ActiveCovs_store lambda_seq




