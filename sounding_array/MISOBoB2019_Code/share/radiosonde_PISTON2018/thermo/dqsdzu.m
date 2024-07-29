function dqsatdzu = dqsdzu(p,Temp)

%  dqsdzu = dqs/dz [K/m] on a dry adiabat
% p [Pa], Temp [degree C]

  global Rd Cp
  dqsatdzu = ((Rd*(Temp+273.15)/Cp) * dqsdT(p,Temp) + p .* dqsdp(p,Temp))./Hscale(Temp+273.15);
% 2017-05-10           ^SPdeS                                                         ^SPdeS