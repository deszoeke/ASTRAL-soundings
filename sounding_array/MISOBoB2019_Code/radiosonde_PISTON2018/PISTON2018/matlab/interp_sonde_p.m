function [hght, temp, rhum, uwnd, vwnd, shum, thta, thte, dwpt] = interp_sonde_p(sonde, pint)
% [hght, temp, rhum, uwnd, vwnd, shum, thta, thte] = interp_sonde_p(sonde, pint)
% Interpolate sonde quantities to pressure levels pint.
% 
% Simon de Szoeke, PISTON 2018

[up, uorder] = unique(sonde.Pressure);

hght = interp1(up, sonde.Height(uorder),            pint);
temp = interp1(up, sonde.Temperature(uorder),       pint);
dwpt = interp1(up, sonde.Dew_Point(uorder),         pint);
rhum = interp1(up, sonde.Relative_Humidity(uorder), pint);

shum = interp1(up, sonde.Specific_Humidity(uorder), pint);
thta = interp1(up, sonde.Potential_Temp(uorder),    pint);
thte = interp1(up, sonde.Equiv_Pot_Temp(uorder),    pint);

uwnd = interp1(up, sonde.u(uorder),                 pint);
vwnd = interp1(up, sonde.v(uorder),                 pint);

end
