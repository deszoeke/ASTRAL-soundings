function Tl = Tlcl(p,Temp,qv)
% Tl=Tlcl(p[Pa], Temp[K], qv[kg/kg])
%  Temperature at the LCL
%   From Bolton, 1980, MWR, 108, 1046-1053.

Rd = 287.04;
Rv = 461.5;

ev = p.*qv./(Rd/Rv + qv);
Tl = 2840./(3.5*log(Temp) - log(.01*ev) - 4.805) + 55;
