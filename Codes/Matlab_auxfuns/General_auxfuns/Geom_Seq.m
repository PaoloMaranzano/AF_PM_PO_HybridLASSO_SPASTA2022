function [geom_seq] = Geom_Seq(min,max,n,zero_flag)

% if isempty(base)
%     geom_seq = exp(linspace(log(min),log(max),n));
% else
%     log_min_base = log(min) / log(base);
%     log_max_base = log(max) / log(base);
%     geom_seq = base.^(linspace(log_min_base,log_max_base,n));
% end

geom_seq = exp(linspace(log(min),log(max),n));


if zero_flag == 1
    geom_seq = [0,geom_seq];
end

return
end