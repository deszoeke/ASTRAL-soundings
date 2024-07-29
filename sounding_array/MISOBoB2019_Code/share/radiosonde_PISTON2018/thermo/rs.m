function rsat = rs(p,T)

%  rs(p,T) is saturation mixing ratio based on Wexler's formula for es
%   with enhancement factor (see es.m).
%   p [Pa], T [degrees C], rs [kg/kg]
%   From A. L. Buck 1981: JAM, 20, 1527-1532.
%   SPdeS 7 July 2004

  global Rd Rv
  esat = es(T,p);
  rsat = (Rd/Rv)*esat./(p-esat);
