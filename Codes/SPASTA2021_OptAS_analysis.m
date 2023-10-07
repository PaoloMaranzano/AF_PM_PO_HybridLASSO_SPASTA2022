%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Dense grid %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
clc

%% %%%%% Loading results/estimates
addpath(genpath('../../../../SPASTA2021/'));
if 0
    %%% r0 results: may 2021
    load('Output/Output_11052021/CV_part.mat')
    load('Output/Output_11052021/Full_model_NonStd.mat')
    load('Output/Output_11052021/LASSO.mat')
    load('Output/Output_11052021/randK_CV.mat')
else
    %%% r1 results: october 2021
    load('Output/Output_0810_dense_grid/CV_part.mat')
    load('Output/Output_0810_dense_grid/output_PMLE_stream.mat')
    load('Output/Output_0810_dense_grid/LASSO.mat')
    load('Output/Output_0810_dense_grid/randK_CV.mat')
end


if 0
    %%%% LASSO output
    [LASSO_RMSE] = LASSO_RMSE_extraction(FitInfo_LASSO,{'NO_2','PM_{10}','PM_{2.5}'},...
        1,0,1,0);
    close all;
    %%%% Mapping LASSO to AS
    [Map_LASSO] = Map_From_Lambda_To_AS(beta_LASSO,Lambda_seq);
    %%%% Optimal AS
    optAS_RMSE_Rand = DSTEM_CV_optAS(randK_CV_metrics,'Random K-fold CV','RMSE',...
        cell2mat(obj_stem_data.stem_varset_p.Y_stds),...
        {'NO_{2}','PM_{10}','PM_{2.5}'},0,0,Map_LASSO,0);
    close all;
    optAS_MAE_Rand = DSTEM_CV_optAS(randK_CV_metrics,'Random K-fold CV','MAE',...
        cell2mat(obj_stem_data.stem_varset_p.Y_stds),...
        {'NO_{2}','PM_{10}','PM_{2.5}'},1,1,Map_LASSO,1);
    close all;
end



if 0
    %% %%%%% Optimal model estimation
    optAS_idx = 1;
    clearvars Setup
    Setup.log_transform = 0;
    Setup.standardize = 1;
    Setup.VarCov_compute = 1;
    Setup.LogLik_compute = 1;
    M1_estim = DSTEM_HDGM_reduced_estim(obj_stem_model_full,Setup,ActiveSets,optAS_idx);
    save('M1_1205.mat','M1_estim')
end

if 1
    %% %%%%% Optimal model estimation
    optAS_idx = 34;
    clearvars Setup
    Setup.log_transform = 0;
    Setup.standardize = 1;
    Setup.VarCov_compute = 1;
    Setup.LogLik_compute = 1;
    M34_estim = DSTEM_HDGM_reduced_estim(obj_stem_model_full,Setup,ActiveSets,optAS_idx);
    save('M34_1305.mat','M34_estim')
end

if 1
    %% %%%%% Optimal model estimation
    optAS_idx = 44;
    clearvars Setup
    Setup.log_transform = 0;
    Setup.standardize = 1;
    Setup.VarCov_compute = 1;
    Setup.LogLik_compute = 1;
    M44_estim = DSTEM_HDGM_reduced_estim(obj_stem_model_full,Setup,ActiveSets,optAS_idx);
    save('M44_1305.mat','M44_estim')
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Large grid %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('Output/Output_09052021/CV_part.mat')
load('Output/Output_09052021/Full_model.mat')
load('Output/Output_09052021/LASSO.mat')
load('Output/Output_09052021/randK_CV.mat')

if 0
    %% %%%%% Optimal model estimation
    optAS_idx = 1;
    clearvars Setup
    Setup.log_transform = 0;
    Setup.standardize = 1;
    Setup.VarCov_compute = 1;
    Setup.LogLik_compute = 1;
    M0_estim = DSTEM_HDGM_reduced_estim(obj_stem_model_full,Setup,ActiveSets,optAS_idx);
    save('M0_1205.mat','M0_estim')
end










clear
clc

load('Output/Output_11052021/LASSO.mat')
load('Output/Output_11052021/randK_CV.mat')
load('Output/Output_12052021/M1_1205.mat')
load('Output/Output_12052021/M34_1205.mat')
load('Output/Output_12052021/M44_1205.mat')

[M1_MLE,M1_perf_met,M1_setup] = DSTEM_extraction(M1_estim,[],[]);
[M34_MLE,M34_perf_met,M34_setup] = DSTEM_extraction(M34_estim,[],[]);
[M44_MLE,M44_perf_met,M44_setup] = DSTEM_extraction(M44_estim,[],[]);

[M1_CV_metrics] = DSTEM_ASi_metric_extraction(M1_estim,randK_CV_metrics,1)
[M34_CV_metrics] = DSTEM_ASi_metric_extraction(M34_estim,randK_CV_metrics,34)
[M44_CV_metrics] = DSTEM_ASi_metric_extraction(M44_estim,randK_CV_metrics,44)


[LASSO_RMSE] = LASSO_RMSE_extraction(FitInfo_LASSO,{'NO_2','PM_{10}','PM_{2.5}'},0,0,0,0);
[Map_LASSO] = Map_From_Lambda_To_AS(beta_LASSO,Lambda_seq);
Lambda_sel = Map_LASSO.Tab_AS([44,34,1],:).min_Lambda;
for m = 1:length(Lambda_sel)
    LASSO_metrics(m,:) = LASSO_RMSE.Trace.Trace_tab(find(LASSO_RMSE.Trace.Trace_tab.Lambda == Lambda_sel(m)),:)
end

load('Output/Output_09052021/LASSO.mat')
load('Output/Output_09052021/randK_CV.mat')
load('Output/Output_12052021/M0_1205.mat')
[M0_MLE,M0_perf_met,M0_setup] = DSTEM_extraction(M0_estim,[],[]);
[M0_CV_metrics] = DSTEM_ASi_metric_extraction(M0_estim,randK_CV_metrics,1)

CV_metrics = array2table([M44_CV_metrics.Avg_RMSE_std , M44_CV_metrics.Avg_MAE_std , M44_CV_metrics.Avg_R2;
    M34_CV_metrics.Avg_RMSE_std , M34_CV_metrics.Avg_MAE_std , M34_CV_metrics.Avg_R2;
    M1_CV_metrics.Avg_RMSE_std , M1_CV_metrics.Avg_MAE_std , M1_CV_metrics.Avg_R2;
    M0_CV_metrics.Avg_RMSE_std , M0_CV_metrics.Avg_MAE_std , M0_CV_metrics.Avg_R2;])
CV_metrics.Properties.VariableNames = {'Avg_RMSE','Avg_MAE','Avg_R2'};
CV_metrics.Model = {'M44','M34','M1','M0'}';

[LASSO_RMSE] = LASSO_RMSE_extraction(FitInfo_LASSO,{'NO_2','PM_{10}','PM_{2.5}'},0,0,0,0);
[Map_LASSO] = Map_From_Lambda_To_AS(beta_LASSO,Lambda_seq);
Lambda_sel = Map_LASSO.Tab_AS([1],:).min_Lambda;
LASSO_metrics(4,:) = LASSO_RMSE.Trace.Trace_tab(find(LASSO_RMSE.Trace.Trace_tab.Lambda == Lambda_sel(1)),:)

LASSO_metrics.Model = {'M44','M34','M1','M0'}';

if 0
    save('Final_models_metrics.mat','M0_CV_metrics','M1_CV_metrics','M34_CV_metrics','M44_CV_metrics','LASSO_metrics','CV_metrics')
end

if 0
    load('Final_models_metrics.mat')
    CV_metrics
    LASSO_metrics
end



t1 = outerjoin(M44_MLE.Reg_pars.Beta_destd_tab.NO2_ground,...
    M34_MLE.Reg_pars.Beta_destd_tab.NO2_ground,'Keys','Var','type','left','MergeKeys',1);
t2 = outerjoin(t1,...
    M1_MLE.Reg_pars.Beta_destd_tab.NO2_ground,'Keys','Var','type','left','MergeKeys',1);
coef_NO2 = outerjoin(t2,...
    M0_MLE.Reg_pars.Beta_destd_tab.NO2_ground,'Keys','Var','type','left','MergeKeys',1);
coef_NO2.Properties.VariableNames = {'Variable',...
    'Coef_M44','SE_coef_M44','abst_coef_M44','Sign_coef_M44',...
        'Coef_M34','SE_coef_M34','abst_coef_M34','Sign_coef_M34',...
            'Coef_M1','SE_coef_M1','abst_coef_M1','Sign_coef_M1',...
                'Coef_M0','SE_coef_M0','abst_coef_M0','Sign_coef_M0'};
            

t1 = outerjoin(M44_MLE.Reg_pars.Beta_destd_tab.PM10_ground,...
    M34_MLE.Reg_pars.Beta_destd_tab.PM10_ground,'Keys','Var','type','left','MergeKeys',1);
t2 = outerjoin(t1,...
    M1_MLE.Reg_pars.Beta_destd_tab.PM10_ground,'Keys','Var','type','left','MergeKeys',1);
coef_PM10 = outerjoin(t2,...
    M0_MLE.Reg_pars.Beta_destd_tab.PM10_ground,'Keys','Var','type','left','MergeKeys',1);
coef_PM10.Properties.VariableNames = {'Variable',...
    'Coef_M44','SE_coef_M44','abst_coef_M44','Sign_coef_M44',...
    'Coef_M34','SE_coef_M34','abst_coef_M34','Sign_coef_M34',...
    'Coef_M1','SE_coef_M1','abst_coef_M1','Sign_coef_M1',...
    'Coef_M0','SE_coef_M0','abst_coef_M0','Sign_coef_M0'};


t1 = outerjoin(M44_MLE.Reg_pars.Beta_destd_tab.PM2_5_ground,...
    M34_MLE.Reg_pars.Beta_destd_tab.PM2_5_ground,'Keys','Var','type','left','MergeKeys',1);
t2 = outerjoin(t1,...
    M1_MLE.Reg_pars.Beta_destd_tab.PM2_5_ground,'Keys','Var','type','left','MergeKeys',1);
coef_PM25 = outerjoin(t2,...
    M0_MLE.Reg_pars.Beta_destd_tab.PM2_5_ground,'Keys','Var','type','left','MergeKeys',1);
coef_PM25.Properties.VariableNames = {'Variable',...
    'Coef_M44','SE_coef_M44','abst_coef_M44','Sign_coef_M44',...
    'Coef_M34','SE_coef_M34','abst_coef_M34','Sign_coef_M34',...
    'Coef_M1','SE_coef_M1','abst_coef_M1','Sign_coef_M1',...
    'Coef_M0','SE_coef_M0','abst_coef_M0','Sign_coef_M0'};


writetable(coef_NO2,'NO2_models.xlsx')
writetable(coef_PM10,'PM10_models.xlsx')
writetable(coef_PM25,'PM25_models.xlsx')





%%%%%%%%% Aggregate coefficients
clear
clc
load('Output/Output_12052021/M34_1205.mat')
[M34_MLE,M34_perf_met,M34_setup] = DSTEM_extraction(M34_estim,[],[]);
names = M34_MLE.Reg_pars.Reg_names_tab.Variable;
idx_names = startsWith(string(names),'lock1');
beta = [M34_MLE.Reg_pars.Beta_destd_tab.NO2_ground.Coef' , ...
    M34_MLE.Reg_pars.Beta_destd_tab.PM10_ground.Coef', ...
    M34_MLE.Reg_pars.Beta_destd_tab.PM2_5_ground.Coef']';
beta_lock = beta(idx_names);
beta_lock_names = names(idx_names);
beta_lock_NO2_vc = M34_MLE.Reg_pars.VarCov_reg_destd.NO2_ground(startsWith(string(M34_MLE.Reg_pars.Reg_names{1}),'lock1'),...
    startsWith(string(M34_MLE.Reg_pars.Reg_names{1}),'lock1'));
coef_NO2_lock = coef_NO2(startsWith(string(coef_NO2.Variable),'lock1'), ...
    {'Variable','Coef_M34','SE_coef_M34'});

poll = {'NO2_ground','PM10_ground','PM2_5_ground'};
for i = 1:3
    names_poll{i} = M34_MLE.Reg_pars.Reg_names{1,i};
    idx_names{i} = startsWith(string(names_poll{i}),'lock1');
    beta_lock{i} = M34_MLE.Reg_pars.Beta_destd_tab.(poll{i})(idx_names{i},1:2);
    beta_lock_vc{i} = M34_MLE.Reg_pars.VarCov_reg_destd.(poll{i})(idx_names{i},idx_names{i});
end

for i = 1:3
    idx_vars_aggr = [1 , 2 , 5];
    c_traf_metr = sum(table2array(beta_lock{i}(idx_vars_aggr,'Coef')));
    se_traf_metr = sqrt(sum(table2array(beta_lock_vc{i}(idx_vars_aggr,idx_vars_aggr)),'all'));
    pv_traf_metr = 1 - normcdf(abs(c_traf_metr/se_traf_metr));
    mt = array2table([c_traf_metr , se_traf_metr, pv_traf_metr]);
    idx_vars_aggr = [1 , 2];
    c_traf = sum(table2array(beta_lock{i}(idx_vars_aggr,'Coef')));
    se_traf = sqrt(sum(table2array(beta_lock_vc{i}(idx_vars_aggr,idx_vars_aggr)),'all'));
    pv_traf = 1 - normcdf(abs(c_traf/se_traf));
    tr = array2table([c_traf , se_traf, pv_traf]);
    idx_vars_aggr = [1 , 4];
    c_rur = sum(table2array(beta_lock{i}(idx_vars_aggr,'Coef')));
    se_rur = sqrt(sum(table2array(beta_lock_vc{i}(idx_vars_aggr,idx_vars_aggr)),'all'));
    pv_rur = 1 - normcdf(abs(c_rur/se_rur));
    rur = array2table([c_rur , se_rur, pv_rur]);
    aggr_effects = [mt ; tr ; rur ];
    aggr_effects.Properties.RowNames = {'Metropolitan traffic','Extra-urban traffic',...
        'Rural background'};
    aggr_effects.Properties.VariableNames = {'Coefficent','SE','PV'};
    aggr_coef{i} = aggr_effects;
end

 




if 0
    
    print(Red_model_estim)
    print(obj_stem_model_full)
    
    
    
    DSTEM_residuals_analysis(obj_stem_model_full)
    
    
    writetable(Full_model_MLE.Reg_pars.Beta_destd_tab.NO2_ground,'NO2_full.xlsx')
    writetable(Full_model_MLE.Reg_pars.Beta_destd_tab.PM10_ground,'PM10_full.xlsx')
    writetable(Full_model_MLE.Reg_pars.Beta_destd_tab.PM2_5_ground,'PM25_full.xlsx')
    writetable(Red_model_MLE.Reg_pars.Beta_destd_tab.NO2_ground,'NO2_red.xlsx')
    writetable(Red_model_MLE.Reg_pars.Beta_destd_tab.PM10_ground,'PM10_red.xlsx')
    writetable(Red_model_MLE.Reg_pars.Beta_destd_tab.PM2_5_ground,'PM25_red.xlsx')
    
    
    
    
    [Full_model_MLE.SpatTime_pars.Theta_Z Model_1_MLE.SpatTime_pars.Theta_Z]
    
    
    avg_RMSE = mean(cell2mat(randK_CV_metrics.RMSE_ASi{1,1}) + ...
        cell2mat(randK_CV_metrics.RMSE_ASi{1,2}) + ...
        cell2mat(randK_CV_metrics.RMSE_ASi{1,3}))
    
    std(avg_RMSE,1)
    
end






%%%%%% Timing analysis
time_end - time_begin
(time_RandKfoldCV_end - time_RandKfoldCV_begin)/10


