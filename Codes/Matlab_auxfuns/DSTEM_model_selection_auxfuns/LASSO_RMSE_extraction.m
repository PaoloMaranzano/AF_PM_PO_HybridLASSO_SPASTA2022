%% %%%%%%%%%% LASSO_RMSE_extraction
%%% Extract the RMSE of each active set in FitInfo_LASSO (DSTEM_SepLASSO).

function [LASSO_RMSE] = LASSO_RMSE_extraction(FitInfo_LASSO,nameY,...
    plot_trace,save_plot_trace,plot_Y,save_plot_Y)

if 0
    nameY = {'NO_2','PM_{10}','PM_{2.5}'};
    plot_trace = 1;
    save_plot_trace = 0;
    plot_Y = 1;
    save_plot_Y = 0;
end



nY = length(FitInfo_LASSO{1,1});
for p = 1:length(FitInfo_LASSO{1,1})
    for i = 1:length(FitInfo_LASSO)
        lam(i) = FitInfo_LASSO{1,i}{1,p}.Lambda;
        lam_RMSE(i) = sqrt(FitInfo_LASSO{1,i}{1,p}.MSE);
    end
    tab_lambda = array2table([lam' , lam_RMSE']);
    tab_lambda.Properties.VariableNames = {'Lambda','RMSE'};
    LASSO_RMSE.Y{1,p}.Tab = tab_lambda;
    
    [val,pos] = min(LASSO_RMSE.Y{1,p}.Tab.RMSE);
    LASSO_RMSE.Y{1,p}.Pos_Lambda_min_RMSE = pos;
    LASSO_RMSE.Y{1,p}.Lambda_min_RMSE = LASSO_RMSE.Y{1,p}.Tab.Lambda(pos);
    LASSO_RMSE.Y{1,p}.Min_RMSE = val;
end

for p = 1:length(FitInfo_LASSO{1,1})
    RMSE_p(:,p) = LASSO_RMSE.Y{1,p}.Tab.RMSE;
end
Trace_RMSE = nanmean(RMSE_p,2);
LASSO_RMSE.Trace.Trace_tab = array2table([LASSO_RMSE.Y{1,1}.Tab.Lambda , Trace_RMSE]);
LASSO_RMSE.Trace.Trace_tab.Properties.VariableNames = {'Lambda','Avg_RMSE'};
[v,p] = min(LASSO_RMSE.Trace.Trace_tab.Avg_RMSE);
LASSO_RMSE.Trace.Min_RMSE = v; 
LASSO_RMSE.Trace.Pos_Lambda_min_RMSE = p;
LASSO_RMSE.Trace.Lambda_min_MSE = LASSO_RMSE.Trace.Trace_tab.Lambda(p); 

if plot_trace == 1
    plot(LASSO_RMSE.Trace.Trace_tab.Lambda,LASSO_RMSE.Trace.Trace_tab.Avg_RMSE);
    xline(LASSO_RMSE.Trace.Lambda_min_MSE,'r');
    yline(LASSO_RMSE.Trace.Min_RMSE,'r');
    xlabel('Lambda')
    ylabel('RMSE')
    title('Average RMSE from separated LASSO')
    if save_plot_trace == 1
        saveas(gcf,'Avg_RMSE_LASSO.png')
    end
end

if plot_Y == 1
    figure()
    for p = 1:nY
        subplot(1,nY,p);
        plot(LASSO_RMSE.Y{1,p}.Tab.Lambda,LASSO_RMSE.Y{1,p}.Tab.RMSE);
        xline(LASSO_RMSE.Y{1,p}.Lambda_min_RMSE,'r');
        yline(LASSO_RMSE.Y{1,p}.Min_RMSE,'r');
        xlabel('Lambda')
        ylabel('RMSE')
        title(['sepLASSO RMSE for ' nameY{p}])
    end
end

return;
end