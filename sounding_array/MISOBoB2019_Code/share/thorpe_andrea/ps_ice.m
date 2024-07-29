function esi = ps_ice(T)
% Saturation vapor pressure over ice given input temperature in Celcius
% Output is in Pascals
    esi = 610.8*exp((21.875 .* T) ./ (T + 265.5));
end

% % Values from a table in a paper
% T = [0.01, 0, -10:-10:-100];
% pi = [611.657, 611.154, 259.923, 103.276, 38.0239, 12.8486, ...
%       3.94017, 1.08204, 0.261893, 0.0548068, 0.00968832, ...
%       0.00140580];
% 
% % Now test the formula from 
% % http://biomet.ucdavis.edu/conversions/HumCon.htm
% tgrid = -100:0.1:0;
% %% plot
% figure; hold on
% scatter(T,pi)
% plot(tgrid,esi)
% 
% % THIS WORKS. 

