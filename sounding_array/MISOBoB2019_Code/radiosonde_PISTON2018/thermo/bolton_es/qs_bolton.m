function qsat = qs_bolton(p,Temp)

%  qs_bolton(p,Temp) is saturation *mixing ratio* based on Wexler's formula for es.
%   From Bolton, 1980, MWR, 108, 1046-1053.
  
  global Rd Rv
  esat = es_bolton(Temp);
  qsat = (Rd/Rv)*esat./(p-esat);
