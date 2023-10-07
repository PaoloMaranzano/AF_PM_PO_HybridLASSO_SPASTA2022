function [Model_CV_metrics] = DSTEM_ASi_metric_extraction(Model_estim,CV_metrics,AS_idx)

if 0
    AS_idx = 1;
    Model_estim = M0_estim;
end

nY = Model_estim.stem_data.stem_varset_p.nvar;
sigmaY = cell2mat(Model_estim.stem_data.stem_varset_p.Y_stds);

for p = 1:nY
    RMSE_ASi = cell2mat(CV_metrics.RMSE_ASi{1,p});
    RMSE(:,p) = RMSE_ASi(:,AS_idx);
    RMSE_std(:,p) = RMSE(:,p) / sigmaY(p);
    MAE_ASi = cell2mat(CV_metrics.MAE_ASi{1,p});
    MAE(:,p) = MAE_ASi(:,AS_idx);
    MAE_std(:,p) = MAE(:,p) / sigmaY(p);
    R2_ASi = cell2mat(CV_metrics.R2_ASi{1,p});
    R2(:,p) = R2_ASi(:,AS_idx);
end

Model_CV_metrics.Avg_RMSE_std = mean(RMSE_std,'all');
Model_CV_metrics.Avg_RMSE = mean(RMSE,'all');
Model_CV_metrics.Avg_MAE_std = mean(MAE_std,'all');
Model_CV_metrics.Avg_MAE = mean(MAE,'all');
Model_CV_metrics.Avg_R2 = mean(R2,'all');

end