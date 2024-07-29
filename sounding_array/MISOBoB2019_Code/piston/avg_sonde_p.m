function [hght, temp, dwpt, rhum, shum, thta, thte, uwnd, vwnd] = avg_sonde_p(sonde, pint)
% [hght, temp, rhum, uwnd, vwnd, shum, thta, thte] = avg_sonde_p(sonde, pint)
% Bin-average sonde quantities in structure sonde between pressure levels pint.
% 
% Simon de Szoeke, PISTON 2018

[up, uorder] = unique(sonde.Pressure);

hght = binavg(pint, up, sonde.Height(uorder)            );
temp = binavg(pint, up, sonde.Temperature(uorder)       );
dwpt = binavg(pint, up, sonde.Dew_Point(uorder)         );
rhum = binavg(pint, up, sonde.Relative_Humidity(uorder) );

shum = binavg(pint, up, sonde.Specific_Humidity(uorder) );
thta = binavg(pint, up, sonde.Potential_Temp(uorder)    );
thte = binavg(pint, up, sonde.Equiv_Pot_Temp(uorder)    );

uwnd = binavg(pint, up, sonde.u(uorder)                 );
vwnd = binavg(pint, up, sonde.v(uorder)                 );
