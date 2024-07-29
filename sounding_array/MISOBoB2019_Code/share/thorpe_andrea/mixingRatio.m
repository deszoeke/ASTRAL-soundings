function [r] = mixingRatio(q)
% Get mixing ratio from specific humidity. Output units will be the same as
% the input units

    r = q ./ (1 - q);

end

