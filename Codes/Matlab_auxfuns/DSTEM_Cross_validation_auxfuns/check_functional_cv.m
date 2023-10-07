function [Y_long_vals,Y_mat_diff] = check_functional_cv(CV_part,data_long,data_table)

for fold = 1:CV_part.K
    mat(:,fold) = CV_part.Y_long_cv{1,fold}{1,1};
end
Y_long_vals = [data_long.Y , mat];

for fold = 1:CV_part.K
    Y_split = splitapply( @(x){x}, Y_long_vals(:,fold+1),data_long.gId );
    Y_fun = cellfun(@transpose,Y_split,'UniformOutput',false);
    for row = 1:size(data_table.Y,1)
        mat2(row,:) = Y_fun{row,1} - data_table.Y{row,1};
    end
    Y_mat_diff{1,fold} = mat2;
end

end