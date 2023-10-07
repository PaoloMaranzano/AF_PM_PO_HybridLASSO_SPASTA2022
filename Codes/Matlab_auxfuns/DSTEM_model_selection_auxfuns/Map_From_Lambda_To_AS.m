
function [Map] = Map_From_Lambda_To_AS(beta_LASSO,Lambda_seq)

beta_LASSO_binary = beta_LASSO ~= 0;
AS = unique(beta_LASSO_binary', 'rows')';
Lambda_AS = nan(1,size(beta_LASSO_binary,2));
for j = 1:size(AS,2)
    Lambda_AS_pos = mean(beta_LASSO_binary == AS(:,j),1) == 1;
    Lambda_AS(Lambda_AS_pos) = j;
end

tab_lambda = array2table([Lambda_seq' , Lambda_AS']);
tab_lambda.Properties.VariableNames = {'Lambda','AS_idx'};

tab_AS = grpstats(tab_lambda,'AS_idx',{'min','mean','median','max'},'DataVars','Lambda');

Map.Tab_lambda = tab_lambda;
Map.Tab_AS = tab_AS;

    return;
end