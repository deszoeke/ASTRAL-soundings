function gamma=moistad(T,p)
% Gamma_m=moistad(T,p)
% Moist adaibatic lapse rate of _potential_ temperature [K/m]
% T in degrees C, p in Pa
% from Wood and Bretherton 2006 (also Rogers and Yau, but with mixing ratio
% instead of specific humidity).
%
% Simon de Szoeke, 2011

global g Cp L Rv Rd

qsat=qs(p,T);
gamma=(g/Cp).*(1-   ( 1+L.*qsat./(Rd.*(T+273.15)) ) ...
                  ./( 1+( (L./(T+273.15)).^2 .*qsat./(Cp*Rv) ) )     );
return