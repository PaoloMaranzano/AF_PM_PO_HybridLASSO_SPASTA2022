%% %%%%% Stratified K-fold cross-validation for spatio-temporal models
if ~exist('Ground')
    inter_lock_type = 1;
    run('..\..\..\Data\DSTEM_data\Reshape_to_HDGM.m');
    Ground.poll = {'NO2','PM10','PM2_5'};
end
clearvars -except crossval_step data Ground log_transform poll standardization

save_cv = 0;

debug_cv = 0;
if debug_cv == 1
    run = 1;
    p = 1;
    Stratum = 1;
end

k = 8;

%%% Partitioning data using a Stratified K-fold CV algorithm
cv_part = CVpart_strat_Kfold(Ground,'Tipology_rec',k);
cv_part = CVpart_strat_Kfold_CondStat(Ground,'Tipology_rec',k);
cv_part = CVpart_random_Kfold(Ground,k);
cv_part = CVpart_strat_boot_CondStat(Ground,'Tipology_rec',k,0.20);
cv_part = CVpart_LeaveOneStatOut(Ground);
%%% Check partitioning
p = 1
tabulate(cv_part.mat_fold{1,p}(:))
%%% Check overlapping
[a,b,c,d] = check_overlap_cv(cv_part,'NO2');
%%% Check proportions
check_prop_strf_cv(cv_part,'PM2_5',8,'Tipology_rec')



