%% %%%%%%%%%% LASSO_MSE_extraction
%%% Extract the MSE of each active set in FitInfo_LASSO (DSTEM_SepLASSO).

function [LASSO_RMSE] = LASSO_RMSE_extraction(FitInfo_LASSO,plot_trace,save_plot_trace)

for p = 1:length(FitInfo_LASSO{1,1})
    for i = 1:length(FitInfo_LASSO)
        lam(i) = FitInfo_LASSO{1,i}{1,p}.Lambda;
        lam_RMSE(i) = sqrt(FitInfo_LASSO{1,i}{1,p}.MSE);
    end
    tab_lambda = array2table([lam' , lam_RMSE']);
    tab_lambda.Properties.VariableNames = {'Lambda','RMSE'};
    LASSO_RMSE{1,p}.Tab = tab_lambda;
    
    [val,pos] = min(LASSO_RMSE{1,p}.Tab.RMSE);
    LASSO_RMSE{1,p}.Pos_Lambda_min_RMSE = pos;
    LASSO_RMSE{1,p}.Lambda_min_RMSE = LASSO_RMSE{1,p}.Tab.Lambda(pos);
    LASSO_RMSE{1,p}.Min_RMSE = val;
end

for p = 1:length(FitInfo_LASSO{1,1})
    RMSE_p(:,p) = LASSO_RMSE{1,p}.Tab.RMSE;
end
Trace_RMSE = sum(RMSE_p,2);
temp_t = array2table([LASSO_RMSE{1,1}.Tab.Lambda , Trace_RMSE]);
LASSO_Trace_RMSE.Trace.Trace_tab = temp_t;
LASSO_RMSE.Trace.Trace_tab.Properties.VariableNames = {'Lambda','Trace_RMSE'};
[v,p] = min(LASSO_Trace_RMSE.Trace.Trace_tab.Trace_RMSE);
LASSO_Trace_RMSE.Trace.Min_RMSE = v; 
LASSO_Trace_RMSE.Trace.Pos_Lambda_min_RMSE = p;
LASSO_RMSE.Trace.Lambda_min_MSE = LASSO_RMSE.Trace.Trace_tab.Trace_RMSE(p); 

if plot_trace == 1
    plot(LASSO_RMSE.Trace_tab.Lambda,LASSO_RMSE.Trace_tab.RMSE);
    xline(LASSO_RMSE.Trace.Lambda_min_MSE,'r');
    yline(LASSO_Trace_RMSE.Trace.Min_RMSE,'r');
    xlabel('Lambda')
    ylabel('RMSE')
    title('Trace RMSE for LASSO')
    if save_plot_trace == 1
        saveas(gcf,'Trace_RMSE_LASSO.png')
    end
end

return;
end