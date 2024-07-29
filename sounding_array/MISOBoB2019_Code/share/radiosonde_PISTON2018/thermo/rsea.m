function rsat = rsea(p,T)

%  rsea(p,T) is saturation mixing ratio over sea water based on an
%   approximation of Wexler's formula for es with enhancement factor
%   (see es.m), and accounting for decreased vapor pressure over saline
%   water.
%   p [Pa], T [degrees C], qsea [kg/kg]
%   From A. L. Buck 1981: JAM, 20, 1527-1532.
%   SPdeS 7. July, 2004
  
  global Rd Rv
  esat = es(T,p) * 0.98; % -2% for salinity (Kraus 1972)
  rsat = (Rd/Rv)*esat./(p-esat);
