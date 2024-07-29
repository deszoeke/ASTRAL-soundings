function esat = es_bolton(Temp)

%  es_bolton(Temp) = 611.2*exp(17.67*(Temp-273.15)/(Temp-273.15+29.5))  [Pa]
%   Saturation vapor pressure (approximate Wexler's formula)
%   From Bolton, 1980, MWR, 108, 1046-1053.
%   REVISED 6. July, 2004 SPdeS to fix for T in Kelvin.
%   es in Pa, T in K.
  
  esat = 611.2*exp(17.67*(Temp-273.15)./(Temp-273.15+243.5));
