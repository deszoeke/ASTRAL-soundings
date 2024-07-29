function qsat = qsea(p,T)

%  qsea(p,T) is saturation specific humidity over sea water based
%   on an approximation of Wexler's formula for es with enhancement
%   factor (see es.m) and reduction of vapor pressure over saline
%   water.
%   p [Pa], T [degrees C], qsea [kg/kg]
%   From A. L. Buck 1981: JAM, 20, 1527-1532.
%   SPdeS 7. July, 2004
  
if ~exist('T','var')
    T=p;
    p=1013+zeros(size(T));
end

  global Rd Rv
  esat = es(T,p) * 0.98; % -2% for salinity (Kraus 1972)
  qsat = (Rd/Rv)*esat./(p+(Rd/Rv-1)*esat);
