function Td=Tdew(e,p)
%   Td=Tdew(e,p)
%   e [Pa] vapor pressure, p pressure [Pa]
%   Dewpoint Td [Celsius] from inversion of es formula of
%   A. L. Buck 1981: JAM, 20, 1527-1532.
%   SPdeS 2 August 2010

  if ~exist('p','var')
    P = 1e3; % default 1000 hPa
  else
    P = p*1e-2; % convert to hPa
  end
  E=e*1e-2;
  
  % e=esat(Td)
  a=17.502;
  b=240.97;
  k=log( E./6.1121./(1.0007+3.46e-6*P) ); % depends only on e,p
  Td=b*k./(a-k);