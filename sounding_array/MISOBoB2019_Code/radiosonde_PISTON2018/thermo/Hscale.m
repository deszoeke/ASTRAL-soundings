function H = Hscale(Temp)
% Hscale(Temp) = Rd*Temp/g [m] is pressure scale height
% Temp [K]
  global Rd g
  H = Rd*Temp/g;