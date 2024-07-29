function Lvap = Lv(Temp)
%  Lv(Temp[K]) = 2.501e6 + (Cpv-Cw)*(Temp-273)  [J/kg]
%   Latent heat of vaporization with temperature correction 
%   From Bolton, 1980, MWR, 108, 1046-1053.

Cpv = 1870; % specific heat of water vapor
Cw  = 4190; % specific heat of liquid water
Lvap = 2.501e6 + (Cpv-Cw)*(Temp-273);