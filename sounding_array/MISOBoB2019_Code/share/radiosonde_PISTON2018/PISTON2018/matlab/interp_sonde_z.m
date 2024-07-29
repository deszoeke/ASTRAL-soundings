function [pres, temp, rhum, uwnd, vwnd, shum, thta, thte, dwpt, wspd, wdir] = interp_sonde_z(sonde, zint)
% [pres, temp, rhum, uwnd, vwnd, shum, thta, thte] = interp_sonde_z(sonde, zint)
% Interpolate sonde quantities to vertical height levels zint.
% 
% Simon de Szoeke, PISTON 2018

[uh, uorder] = unique(sonde.Height);

pres = interp1(uh, sonde.Pressure(uorder),          zint);
temp = interp1(uh, sonde.Temperature(uorder),       zint);
dwpt = interp1(uh, sonde.Dew_Point(uorder),         zint);
rhum = interp1(uh, sonde.Relative_Humidity(uorder), zint);

shum = interp1(uh, sonde.Specific_Humidity(uorder), zint);
thta = interp1(uh, sonde.Potential_Temp(uorder),    zint);
thte = interp1(uh, sonde.Equiv_Pot_Temp(uorder),    zint);

uwnd = interp1(uh, sonde.u(uorder),                 zint);
vwnd = interp1(uh, sonde.v(uorder),                 zint);

wspd = interp1(uh, sonde.Wind_Speed(uorder),        zint);
wdir = interp1(uh, double(sonde.Wind_Direction(uorder)),    zint);

end
