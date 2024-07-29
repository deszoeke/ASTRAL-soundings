function B_la=Planck_la(la,T)
% B_la=Planck_la(lambda[m],T[K]) [W/sr/m^3]
% Computes Planck's blackbody spectral radiance as a function of wavelength lambda [m].
% (c) Simon de Szoeke 2016-04-18

% constants
kb=1.381e-23;   % J/K   Boltzmann
% h=6.6261e-34;   % J s   Planck
c=2.998e8;     % m/s   speed of light
hc=1.98645e-25; % J m   Planck*speed of light

B_la=2*hc*c./((la.^5).*(exp(hc./(la*kb*T))-1)); % W/sr/m^3
%B_nu=2*h*nu.*nu.*nu./(c*c.*(exp(h*nu./(kb*T))-1)); % W/sr/m^2 s