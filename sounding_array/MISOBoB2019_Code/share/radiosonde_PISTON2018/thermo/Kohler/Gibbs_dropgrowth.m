% plot Gibbs energy as a function of droplet radius r
% cd ~/Dropbox/ATS511/Week07

e=612;  % Pa
es=611; % Pa
% S=e/es; % or
rhol=1e3; % kg m^-3
Rv=461; % J/K/kg
T=273.15;
sigma=75.6e-3; % J/m^2

r=[6e-8:2e-8:1.5e-6]'; % radius m
clf
plot([0 1.5],[0 0],'k'); hold on
% subsaturation
S=0.999;
g2= -rhol*Rv*T*log(S)*4/3*pi*r.^3;
g1= 4*pi*sigma*r.^2;
G=g1+g2;
plot(r*1e6,G,'r','linewidth',1.5)

% saturation
S=1;
g2= -rhol*Rv*T*log(S)*4/3*pi*r.^3;
g1= 4*pi*sigma*r.^2;
G=g1+g2;
plot(r*1e6,G,'linewidth',1.5)

% supersaturation
S=1.002;
g2= -rhol*Rv*T*log(S)*4/3*pi*r.^3;
g1= 4*pi*sigma*r.^2;
G=g1+g2;
plot(r*1e6,G,'linewidth',1.5,'color',[.0 .6 .2])

set(gca,'fontsize',16)
ylabel('Gibbs energy (J)')
xlabel('droplet radius (10^{-6} m)')

print -depsc Gibbs_dropgrowth1.eps

axis([0.2 1 -1e-13 2.5e-13])
print -depsc Gibbs_dropgrowth2.eps

% plot critical saturation as a function or radius
r=10.^(-8:.1:-4)'; % m
S=exp(2*sigma./(rhol*Rv*T.*r));
h=semilogx(r*1e6,S,'b',[1e-2 1e2],[1 1],'k'); set(h(1),'linewidth',2)
set(gca,'fontsize',16)
axis tight
xlabel('droplet radius (10^{-6} m)')
ylabel('saturation  \it{S}\rm = \it{e_s}\rm(\it{r}\rm)/\it{e_s}\rm(\infty)')
print -depsc crit_supersat.eps
