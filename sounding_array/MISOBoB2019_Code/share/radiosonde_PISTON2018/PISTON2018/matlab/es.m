function esat = es(T,p,forceC)

%  es(T,p) = is saturation vapor pressure based on Wexler's formula,
%   with enhancement factor for moist air rather than water vapor.
%   The enhancement factor requires a pressure.
%   T [degrees C], p [Pa] (note the reversed input order), es [Pa]
%   From A. L. Buck 1981: JAM, 20, 1527-1532.
%   SPdeS 7 July 2004

  if ~exist('p','var')
    P = 1e3; % default 1000 hPa
  else  
    P = p*1e-2; % convert to hPa
  end
  if ~exist('forceC','var')
    forceC = logical(0);
  end
  if T>200 & ~forceC
    T = T - 273.15;
    warning('Temperature in Kelvin detected.')
  end
  esat = 6.1121*(1.0007 + 3.46e-8*P).*exp((17.502*T)./(240.97 + T));
  esat = esat*1e2; % convert to Pa
