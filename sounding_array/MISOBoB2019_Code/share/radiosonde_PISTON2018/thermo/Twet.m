function Tw=Twet(Tdry,qv,p,niter)
% Tw=Twet(Tdry[C], qv[kg/kg], p[Pa])
% Computes the wet bulb temperature from an isobaric process
% by using heat internal energy -Cp*dT=L*dq to evaporate
% and raise qv to qs(Twb).
%
% Simon de Szoeke 2016-04-07

if ~exist('niter','var')
    niter=5;
end

global Cp Cpv Cw

% initializes Tw, qw
Tw=Tdry;
% problem for temperatures greater than 100 C!
qw=min(qs(p,Tw),qv);

% Steps by 1/4 of saturation deficit and enforces -Cp*dT=L*dq by successive approximations
for iter=1:niter
    dq=0.25*(qs(p,Tw)-qw);
    if dq<0; break; end
    Lvap = 2.501e6 + (Cpv-Cw)*Tw;
    dT=-Lvap/Cp.*dq;
    qw=qw+dq;
    Tw=Tw+dT;
end
