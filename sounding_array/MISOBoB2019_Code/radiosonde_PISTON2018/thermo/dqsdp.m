function dqsatdp = dqsdp(p,Temp)
%  dqsdp(p,Temp) = d(qs)/d(Temp) = -qs(p,Temp)/(p - es(Temp,p))   [Pa^{-1}]
%  p[Pa], Temp[degree C]
  dqsatdp = -qs(p,Temp)./(p - es(Temp,p));
