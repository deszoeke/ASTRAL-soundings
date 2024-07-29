function dTdz = dTdzsat(p , Temp , inflag, outflag)
% DOES NOT WORK!
% (K/m) = dTdzsat(p (Pa), Temp (K), [inflag, [outflag]])
%  Given a temperature (or potential temperature) and pressure, determine
%  dT/dz for a moist adiabat.
%   If inflag is nonzero, Temp input is treated as theta.
%   If outflag is nonzero, dtheta/dz is returned.
  
  thermo_constants
  
  if ~exist('inflag' ,'var'); inflag  = 0; end;
  if ~exist('outflag','var'); outflag = 0; end;
  
  Pi = (p./1e5).^(Rd/Cp);
  if(inflag) % input is theta
    th = Temp
    Temp = th.*Pi
  else       % input is T
    th = Temp./Pi;
  end
  
  dTdz  = -g/Cp./(1 + (Lv(Temp)/Cp)*dqsdT(p,Temp))
%                     ^^^^^^^^^^^Gamma^^^^^^^^^^^
  dTdz = -Lv(Temp)/Cp*dqsdzs(p,Temp) - g/Cp;
  gamma = -dTdz;
  gammadry = g/Cp;
  -(gammadry-gamma)./Pi
  dthdz = (g/Cp - dTdz)./Pi
  
  if(outflag)
    dTdz = dthdz;
  end
  