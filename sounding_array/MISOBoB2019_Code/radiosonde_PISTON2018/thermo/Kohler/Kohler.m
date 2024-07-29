% produce Koehler curve
e=612;  % Pa
es=611; % Pa
% S=e/es; % or
rhol=1e3; % kg m^-3
Rv=461; % J/K/kg
T=273.15;
sigma=75.6e-3; % J/m^2
Mv=18.016;   % molar mass of water

a=2*sigma/(rhol*Rv*T); % size coefficient m

% common solutes
%           molar mass      density
%NaCl       58.44           2.165e3
%(NH4)2 SO4 132.14          1.77e3
%(NH4)2 NO3 80.1            1.73e3
ms=1e-19; % solute mass kg
Ms=132.14;  % molar mass of solute kg/kmol (ballpark)
vH=2; % van't Hoff factor
b=3*vH*ms*Mv/(4*pi*rhol*Ms); % solution coefficient m^3

r=10.^(-8:.02:-5)'; % m
solute=1-b./r.^3;
sizeff=exp(a./r);
plot([1e-8 1e-5],[1 1],'k:'); hold on
hl=plot(r,sizeff,'-b',r,solute,'r--',r,solute.*sizeff);
set(hl(3),'linewidth',1.5,'color',[.05 .5 .1])
set(hl(1:2),'linewidth',0.8)
set(gca,'ylim',[0.985 1.02],'xlim',[2e-8 1e-5],'xscale','log','fontsize',16)
title('Kohler curve (10^{-19} kg ammonium sulfate)')
ylabel('saturation ratio  \it{S}\rm = \it{e}/\it{e_s}\rm(\infty)')
xlabel('droplet radius (m)')
%print -depsc Kohler.eps

% 10^-18 kg for different solute species and masses
clf
%NaCl
ms=1e-18; % solute mass kg
Ms=58.4;  % molar mass of solute kg/kmol
vH=2; % van't Hoff factor
b=3*vH*ms*Mv/(4*pi*rhol*Ms); % solution coefficient m^3
solute=1-b./r.^3;
plot([1e-8 1e-5],[1 1],'k:'); hold on
hl=plot(r,sizeff,'-b',r,solute,'r:',r,solute.*sizeff);
hold on
set(hl(3),'linewidth',1.5,'color',[.05 .5 .5])
set(hl(2),'color',[.05 .5 .5])

%(NH4)2NO3
ms=1e-18; % solute mass kg
Ms=80.1;  % molar mass of solute kg/kmol
vH=2; % van't Hoff factor
b=3*vH*ms*Mv/(4*pi*rhol*Ms); % solution coefficient m^3
solute=1-b./r.^3;
plot([1e-8 1e-5],[1 1],'k:'); hold on
hl(4:5)=plot(r,solute,'r:',r,solute.*sizeff);
hold on
set(hl(5),'linewidth',1.5,'color',[.6 .2 0])
set(hl(4),'color',[.6 .2 0])

%(NH4)2SO4
ms=1e-18; % solute mass kg
Ms=132.14;  % molar mass of solute kg/kmol (ballpark)
vH=2; % van't Hoff factor
b=3*vH*ms*Mv/(4*pi*rhol*Ms); % solution coefficient m^3
solute=1-b./r.^3;
hl(6:7)=plot(r,solute,'r:',r,solute.*sizeff);
set(hl(7),'linewidth',1.5,'color',[.5 .5 0])
set(hl(6),'color',[.5 .5 0])

set(gca,'ylim',[0.99 1.005],'xlim',[1e-7 1e-5],'xscale','log','fontsize',16)
leg=legend(hl([3 5 7]),'10^{-18} kg NaCl','10^{-18} kg (NH_4)_2NO_3','10^{-18} kg (NH_4)_2SO_4');
ylabel('saturation ratio  \it{S}\rm = \it{e}/\it{e_s}\rm(\infty)')
xlabel('droplet radius (m)')
print -depsc Kohler_NaCl_AN_AS.eps

% RY prob 6.3
%NaCl
T=280;
ms=1e-19; % solute mass kg
Ms=58.4;  % molar mass of solute kg/kmol (ballpark)
vH=2; % van't Hoff factor
b1=3*vH*ms*Mv/(4*pi*rhol*Ms); % solution coefficient m^3
solute1=1-b1./r.^3;
S1=solute1.*sizeff;
%(NH4)2SO4
Ms=132.14;  % molar mass of solute kg/kmol (ballpark)
vH=2; % van't Hoff factor
b2=3*vH*ms*Mv/(4*pi*rhol*Ms); % solution coefficient m^3
solute2=1-b2./r.^3;
S2=solute2.*sizeff;
subplot(2,1,1)
hl=plot(r,solute1.*sizeff,r,solute2.*sizeff);
set(hl(:),'linewidth',1.5,'color',[.1 .5 .5])
set(gca,'ylim',[0.985 1.02],'xlim',[2e-8 1e-5],'xscale','log','fontsize',16)
set(hl(2),'color',[.5 .5 0])
legend('10^{-19} kg NaCl','10^{-19} kg (NH_4)_2SO_4',0)
hold on
plot([1e-8 1e-5],[1 1],'k:')
ylabel('saturation ratio  \it{S}\rm = \it{e}/\it{e_s}\rm(\infty)')
xlabel('droplet radius (m)')
subplot(2,1,2)
plot(r,S2-S1)
set(gca,'ylim',[0 1.5],'xlim',[2e-8 1e-5],'xscale','log','fontsize',16)
print -depsc Kohler_RYprob6p3.eps

% exam problem
T=280;
ms=1e-19;
Ms=60;
vH=2; % van't Hoff factor
b=3*vH*ms*Mv/(4*pi*rhol*Ms); % solution coefficient m^3
sigma1=7.5e-2; % J m^2
a1=2*sigma1/(rhol*Rv*T); % size coefficient m
sigma2=3.75e-2;
a2=2*sigma2/(rhol*Rv*T); % size coefficient m
r=10.^(-8:.02:-5)'; % m
curv1=exp(a1./r);
curv2=exp(a2./r);
solute=1-b./r.^3;
semilogx(r,curv1,r,curv2)
hold on
semilogx(r,a1./r+solute,'--',r,a2./r+solute,'--')
plot(r,solute.*curv1,r,solute.*curv2)
plot([1e-8 1e-5],[1 1],'k:')
set(gca,'ylim',[0.985 1.02],'xlim',[2e-8 1e-5],'xscale','log','fontsize',16)
ylabel('saturation ratio  \it{S}\rm = \it{e}/\it{e_s}\rm(\infty)')
xlabel('droplet radius (m)')
