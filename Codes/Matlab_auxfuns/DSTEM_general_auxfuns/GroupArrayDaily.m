function Var = GroupArrayDaily(Varargin,VarName,gID,DatasetName)
    Var = splitapply( @(x){x}, Varargin, gID );
    Var = cellfun(@transpose,Var,'UniformOutput',false);
    Name = DatasetName;
    if i==1 & ~exist(Name)
        output_final2 = table(Var);
        output_final2.Properties.VariableNames{'Var'} = VarName;
    elseif i==1 & exist(Name)
        clear DatasetName
        output_final2 = table(Var);
        output_final2.Properties.VariableNames{'Var'} = VarName;
    else
        output_final2.X_beta_var = Var;
        output_final2.Properties.VariableNames{'X_beta_var'} = VarName;
    end
end