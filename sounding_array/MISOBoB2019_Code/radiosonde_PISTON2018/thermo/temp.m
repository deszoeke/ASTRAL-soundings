function t = temp(p,th)
% Temperature theta*(1e5/p)^(-Rd/Cp)
% p (Pa), theta (K); temperature (K)
  
  global Rd Cp
  t = th.*(1e5./p).^(-Rd/Cp);
