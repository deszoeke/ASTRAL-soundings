function dqsatdzs = dqsdzs(p,Temp)
%  dqsdzs(p,Temp) = dqs/dz [K/m] on a moist adiabat.
%  p [Pa], Temp [C]
%
% dqsatdzs = dqsdzu(p,Temp)./(1 + Gamma(p,Temp));
  dqsatdzs = dqsdzu(p,Temp)./(1 + Gamma(p,Temp));
