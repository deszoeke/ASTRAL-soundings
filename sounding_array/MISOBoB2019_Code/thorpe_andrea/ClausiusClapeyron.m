function [es] = ClausiusClapeyron(T)
% Saturation specific humidity at an input temperature in Kelvin.
% Output is in Pascals. 

% NEVERMIND: Use the August-Roche-Magnus approximation

%     if any(T < 0)
%         disp('Warning: Negative T detected. Check for Kelvin units')
%     end

    L  = 2.5e6;
    Rv = 461;
    
    es = 611 * exp((L/Rv) .* ((1/273) - (1./T)));

%     es = 610.94 * exp((17.625 .* T)./(T + 243.04));
end

