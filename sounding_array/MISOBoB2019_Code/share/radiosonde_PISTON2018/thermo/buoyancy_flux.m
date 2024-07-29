function [bflx, bflx_s]=buoyancy_flux(shf, lhf, T, q, p)
% [bflx, bflx_s]=buoyancy_flux(shf, lhf, T, q, p)  [shf, lhf W/m^2; T °C; q g/kg; p hPa]
% Buoyancy flux from sensible and latent heat flux. Optionally also returns the sensible
% part of the buoyancy flux. bflx = bflx_virt + bflx_s
%
% Simon de Szoeke  2018-09-01

thermo_constants % loads Rd, Rv, Cp, g, etc.

epsi = 0.608; % Rv/Rd - 1
virtfac=1.0+epsi*q/1e3;
Tv = (T+273.15).*virtfac;
rho = p*1e2./(Rd*Tv);
wt = shf./(rho*Cp);
wtv = wt .* virtfac + lhf./(rho.*Lv(T+273.15)) .* (T+273.15)*epsi;
bflx = g*wtv./Tv;
bflx_s = g*wt./Tv;

end