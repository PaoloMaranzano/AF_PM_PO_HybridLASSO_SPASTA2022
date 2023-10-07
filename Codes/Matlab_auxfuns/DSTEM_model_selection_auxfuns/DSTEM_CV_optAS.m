function [optAS_info] = DSTEM_CV_optAS(cv_metrics,cv_name,metric,stdY,nameY,...
    plot_metric,save_p_metric,map_lambda_AS,save_p_map)

if 0
    stdY = cell2mat(obj_stem_data.stem_varset_p.Y_stds);
    cv_metrics = randK_CV_metrics;
    cv_name = 'RandK';
    metric = 'RMSE';
    nameY = {'NO2','PM_{10}','PM_{2.5}'};
    map_lambda_AS = Map_LASSO;
    plot_metric = 0;
end

% Useful quantities
nY = size(cv_metrics.([metric '_ASi']),2);
K = size(cv_metrics.([metric '_ASi']){1,1},1);
if metric == "RMSE"
    CV_metric = cv_metrics.('MSE_ASi');
else
    CV_metric = cv_metrics.([metric '_ASi']);
end

for p = 1:nY
    C = CV_metric{1,p};
    fx=@(x)any(isempty(x));
    ind=cellfun(fx,C);
    C(ind)={nan};
    metric_Y{p} = cell2mat(C);
    if metric == "RMSE"
        metric_Y{p} = sqrt(metric_Y{p});
    end
    if metric ~= "R2"
        metric_Y{p} = metric_Y{p} ./ stdY(p);
    end
    avg_metric(p,:) = mean(metric_Y{p},1);
    if metric ~= "R2"
        [Val,Pos] = min(avg_metric(p,:));
    else
        [Val,Pos] = max(avg_metric(p,:));
    end
    Y_metric_optValue(p) = Val;
    Y_metric_optAS(p) = Pos;
end

metric_Y_3d = reshape(cell2mat(metric_Y),K,[],nY);
trace_metric = nanmean(metric_Y_3d,3);
avg_trace_metric = mean(trace_metric,1);
if metric ~= "R2"
    [Val,Pos] = min(avg_trace_metric);
else
    [Val,Pos] = max(avg_trace_metric);
end
trace_metric_optValue = Val;
trace_metric_optAS = Pos;

if plot_metric == 1
    % Average metric
    figure()
    plot(avg_trace_metric);
    xline(trace_metric_optAS,'r');
    yline(trace_metric_optValue,'r');
    xlabel('Active sets from sparsest (left) to fullest (right)')
    ylabel([metric])
    title(['Average ' metric ' for ' cv_name])
    if save_p_metric == 1
        saveas(gcf,[metric '_' cv_name '.png'])
    end
    % Metric for each response
    figure()
    for p = 1:nY
        if metric ~= "R2"
            [val,pos] = min(avg_metric(p,:));
        else
            [val,pos] = max(avg_metric(p,:));
        end
        subplot(1,nY,p);
        plot(avg_metric(p,:));
        xline(pos,'r');
        yline(val,'r');
        xlabel('Active sets from sparsest (left) to fullest (right)')
        ylabel([metric])
        title(['Average ' metric ' for ' nameY{p}])
    end
    if save_p_metric == 1
        saveas(gcf,['Single_' metric '_' cv_name '.png'])
    end
end

if ~isempty(map_lambda_AS)
    DSTEM_RMSE_lambda = avg_trace_metric(map_lambda_AS.Tab_lambda.AS_idx);
    if metric ~= "R2"
        [v,p] = min(DSTEM_RMSE_lambda);
    else
        [v,p] = max(DSTEM_RMSE_lambda);
    end
    figure()
    plot(map_lambda_AS.Tab_lambda.Lambda,DSTEM_RMSE_lambda);
    xline(map_lambda_AS.Tab_lambda.Lambda(p),'r');
    yline(v,'r');
    xlabel('Lambda')
    ylabel([metric])
    title(['Average ' metric ' for ' cv_name])
    if save_p_map == 1
        saveas(gcf,['Map_lambda_' metric '_' cv_name '.png'])
    end
    
    figure()
    for p = 1:nY
        DSTEM_RMSE_lambda_p = avg_metric(p,map_lambda_AS.Tab_lambda.AS_idx);
        if metric ~= "R2"
            [val,pos] = min(DSTEM_RMSE_lambda_p);
        else
            [val,pos] = max(DSTEM_RMSE_lambda_p);
        end
        subplot(1,nY,p);
        plot(map_lambda_AS.Tab_lambda.Lambda,DSTEM_RMSE_lambda_p);
        xline(map_lambda_AS.Tab_lambda.Lambda(pos),'r');
        yline(val,'r');
        xlabel('Lambda')
        ylabel([metric])
        title(['Average ' metric ' for ' nameY{p}])
    end
    if save_p_map == 1
        saveas(gcf,['Single_map_lambda_' metric '_' cv_name '.png'])
    end
end
    

%%% Output
optAS_info.Metric = metric;
optAS_info.Y_metric_optValue = Y_metric_optValue;
optAS_info.Y_metric_optAS = Y_metric_optAS;
optAS_info.Trace_metric_avg = avg_trace_metric;
optAS_info.Trace_metric_optValue = trace_metric_optValue;
optAS_info.Trace_metric_optAS = trace_metric_optAS;

return;
end