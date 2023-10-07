function [y] = round_n_digits(x,n)
% "Round the matrix x to n decimal places"
if istable(x) == 1
    name_cols = x.Properties.VariableNames;
    name_rows = x.Properties.RowNames;    
    y = table2array(x);
end

y = round(y * 10^n) / 10^n;

if istable(x) == 1 
    y = array2table(y);
    y.Properties.VariableNames = name_cols;
    y.Properties.RowNames = name_rows;
end

end

