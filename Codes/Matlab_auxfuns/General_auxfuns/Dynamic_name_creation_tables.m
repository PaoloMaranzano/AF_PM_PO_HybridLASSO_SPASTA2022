%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Dynamic creation of tables %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%%% Create matrices
A = [1 2 3 4 ;
     5 6 7 8 ;
     9 10 11 12;
     13 14 15 16;]
B = randn(10,4);
%% %%%%% Convert matrices into tables with variable names
T1 = array2table(A);
T1.Properties.VariableNames(1:4) = {'QA','QB','QC','QD'}
T2 = array2table(B);
T2.Properties.VariableNames(1:4) = {'QA','QB','QC','QD'}

%% Dynamically assignment to new tables with sequential names
for i=1:2
genvarname('Q',  num2str(i));
eval(['Q' num2str(i) '= T' num2str(i)])
end