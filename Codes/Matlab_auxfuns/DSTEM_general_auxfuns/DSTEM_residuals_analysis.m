%% %%%%%%%%%% DSTEM_residuals_analysis
%%% Computes and pltos the residuals for an estimated DSTEM model.

function [Residuals] = DSTEM_residuals_analysis(DSTEM_model)

debug = 0;
if debug == 1
	DSTEM_model = obj_stem_model_full;
end

nY = DSTEM_model.stem_data.stem_varset_p.nvar;
nameY = DSTEM_model.stem_data.stem_varset_p.Y_name;
nameY = {'NO_2','PM_{10}','PM_{2.5}'};

Residuals = DSTEM_model.stem_EM_result.res;

figure()
for p = 1:nY
    subplot(nY,1,p)
    plot(Residuals{1,p}')
    title(['Residuals for ' nameY{p}])
    ylabel('Residuals')
    xlabel('Time stamps')
end

return;
end



