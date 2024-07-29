function [pres, temp, rhum, uwnd, vwnd, shum, thta, thte, dwpt, wspd, wdir] = avg_sonde_z_hed(sonde, zint)
% [pres, temp, rhum, uwnd, vwnd, shum, thta, thte, dwpt, wspd, wdir] = avg_sonde_z_hed(sonde, zint)
% Bin-average sonde quantities from structure sonde between vertical height levels zint.
% 
% Simon de Szoeke, PISTON 2018

[uh, uorder] = unique(sonde.Height);

pres = binavg(zint, uh, sonde.Pressure(uorder)          );
temp = binavg(zint, uh, sonde.Temperature(uorder)       );
dwpt = binavg(zint, uh, sonde.Dew_Point(uorder)         );
rhum = binavg(zint, uh, sonde.Relative_Humidity(uorder) );

shum = binavg(zint, uh, sonde.Specific_Humidity(uorder) );
thta = binavg(zint, uh, sonde.Potential_Temp(uorder)    );
thte = binavg(zint, uh, sonde.Equiv_Pot_Temp(uorder)    );

uwnd = binavg(zint, uh, sonde.u(uorder)                 );
vwnd = binavg(zint, uh, sonde.v(uorder)                 );

wspd = binavg(zint, uh, sonde.Wind_Speed(uorder)        );
wdir = binavg(zint, uh, sonde.Wind_Direction(uorder)    );
end
