%% Check overlapping
function [check_cv,diff_by_row,diff_overall,diff_max] = check_overlap_cv(CV_part,Pollutant)

% CV_part = cv_part;
% Pollutant = 'NO2';
poll = CV_part.poll;
p = find(ismember(poll, Pollutant));
K = CV_part.K;
if CV_part.cv_type == "StratBoot"
    % Stratified bootstrap case
    for f = 1:K
        check_cv(:,:,f) = CV_part.mat_fold{1,f}{1,p} == 1 | isnan(CV_part.mat_fold{1,f}{1,p});
    end
else
    % Other cases
    for f = 1:K
        check_cv(:,:,f) = CV_part.mat_fold{1,p} == f | isnan(CV_part.mat_fold{1,p});
    end
end
check_cv = sum(check_cv,3);
diff_max = tabulate(check_cv(check_cv~=K));
diff_by_row = sum((check_cv == K) - isnan(CV_part.Y{1,p}),2);
diff_overall = sum(sum((check_cv == K) - isnan(CV_part.Y{1,p}),2));

return;

end
