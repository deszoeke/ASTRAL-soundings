function qsat = qs(p,T)
%  qs(p,T) is saturation specific humidity based on Wexler's formula for es
%   with enhancement factor (see es.m).
%   p [Pa], T [degrees C], qs [kg/kg]
%   From A. L. Buck 1981: JAM, 20, 1527-1532.
%   SPdeS 7 July 2004

Rd = 287.04;
Rv = 461.5;
RdoRv=Rd/Rv;

esat = es(T,p);
qsat = RdoRv*esat ./ (p + (RdoRv-1)*esat);
