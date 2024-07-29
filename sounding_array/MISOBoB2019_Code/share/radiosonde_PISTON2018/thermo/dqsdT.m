function dqsatdT = dqsdT(p,Temp)

%  dqsdT = d(qs(p, Temp))/d(Temp) = qs(p,Temp)*Lv(Temp)/(Rv*Temp**2)  [K^{-1}]
%   using Wexler and Clausius-Clapeyron formulas
%   p[Pa], T[degrees C]
  
%  273, 273.15 shifts in Temp are deliberate, to offset other
%  inconsistencies in the thermo functions. SPdeS 01.July.2005
  
  global Rv
  dqsatdT = qs(p,Temp).*Lv(Temp+273)./(Rv*(Temp+273.15).^2);
