function gam = Gamma(p,Temp)
%  Gamma(p,Temp) = (Lv(Temp)/Cp)*dqsdT(p,Temp)
%  p [Pa], Temp [degree C]
  global Cp
  gam =  (Lv(Temp+273.)/Cp).*dqsdT(p,Temp);
% 2017-05-10     ^SPdeS, 273. deliberately compensates Lv
