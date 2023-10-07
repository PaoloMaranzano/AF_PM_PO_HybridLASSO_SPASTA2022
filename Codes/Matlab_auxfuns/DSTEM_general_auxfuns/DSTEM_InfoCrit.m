%% %%%%%%%%%% DSTEM_InfoCrit
%%% Computes the ML Information Criteria for an estimated DSTEM model.

function [AIC,AICc,BIC,HQIC] = DSTEM_InfoCrit(DSTEM_model)

%%% Number of parameters
npars = length(DSTEM_model.par_vec());
%%% Number of observations
nobs = size(DSTEM_model.stem_data.Y,1)*size(DSTEM_model.stem_data.Y,2);
%%% Log-likelihood
% LogL = DSTEM_model.stem_EM_result.logL_all(end);
LogL = DSTEM_model.stem_EM_result.logL;
%%% Information criteria
AIC = LogL - 2*npars;
AICc = DSTEM_model.stem_EM_result.AIC + (2*npars^2 + 2*npars)/(nobs - npars - 1);
BIC = LogL - log(nobs)*npars;
HQIC = LogL - 2*npars*log(log(nobs));

end