cd('/Users/sdeszoek/Data/cruises/PISTON_MISOBOB_2018/Soundings/Thompson/matlab')

addpath(genpath('/Users/sdeszoek/Documents/MATLAB/user_tools/stat/'));
addpath(genpath('/Users/sdeszoek/Documents/MATLAB/user_tools/thermo/'));
addpath(genpath('/Users/sdeszoek/Documents/MATLAB/user_tools/graphics/'));
thermo_constants;

update = false;
if update
    % download latest U.Wyo and TGT sounding files:
    rsync_tgt = 'rsync -parv /Volumes/cruiseshare/Soundings/Thompson/ ~/Data/cruises/PISTON_MISOBOB_2018/Soundings/Thompson/';
    fetch_pal = 'curl -s "http://weather.uwyo.edu/cgi-bin/sounding?region=pac&TYPE=TEXT%3ALIST&YEAR=2018&MONTH=09&FROM=0100&TO=3012&STNM=91408" | sed -e "s/<[a-zA-Z\/][^>]*>//g" > ~/Data/cruises/PISTON_MISOBOB_2018/Soundings/Palau/91408_KororPalau/201809.txt &';
    fetch_yap = 'curl -s "http://weather.uwyo.edu/cgi-bin/sounding?region=pac&TYPE=TEXT%3ALIST&YEAR=2018&MONTH=09&FROM=0100&TO=3012&STNM=91413" | sed -e "s/<[a-zA-Z\/][^>]*>//g" > ~/Data/cruises/PISTON_MISOBOB_2018/Soundings/Yap/201809.txt &';
    
    %system(rsync_tgt)
    system(fetch_pal)
    system(fetch_yap)
    % preprocess Yap and Palau with python scripts
    % ../../Yap/read_sondes_and_graphs.py
    % ../../Palau/91408KororPalau/read_sondes_and_graphs.py
    % link Yap_soundings.mat and Palau_soundings.mat to ...Thompson/matlab/

    % bin to consistent pressure coordinates and write to mat file
    read_pbin_all_TGT_sondes;
    read_pbin_all_YapPalau_sondes;
end

plev = [0:25:150 200:50:1000]';
dp = diff(plev);
pmid = plev(1:end-1) + dp/2;
np = length(pmid);

Palau = load('Palau_sondes_pavg.mat');
Palau.lat = 7.34;
Palau.lon = 134.48;

Yap = load('Yap_sondes_pavg.mat');
Yap.lat = 9.5;
Yap.lon = 138.08;

TGT = load('TGT_sondes_pavg.mat');
% ship lat and lon vary

% master 6 h time grid
mday6h = datenum(2018,8,22,12,0,0):0.25:min(now,datenum(2018,10,9));
yday6h = mday6h - datenum(2018,0,0);
nt = length(yday6h);

[T.uwnd, T.vwnd, T.temp, T.thta, T.shum, T.hght] = deal(NaN(np,nt));
T.lat = binavg(yday6h-0.125, TGT.yday, TGT.lat);
T.lon = binavg(yday6h-0.125, TGT.yday, TGT.lon);
for i = 1:np
    % bin average TGT to 6 h
    T.uwnd(i,:) = binavg(yday6h-0.125, TGT.yday, TGT.uwnd(i,:));
    T.vwnd(i,:) = binavg(yday6h-0.125, TGT.yday, TGT.vwnd(i,:));
    T.temp(i,:) = binavg(yday6h-0.125, TGT.yday, TGT.temp(i,:));
    T.thta(i,:) = binavg(yday6h-0.125, TGT.yday, TGT.thta(i,:));
    T.shum(i,:) = binavg(yday6h-0.125, TGT.yday, TGT.shum(i,:));
    T.hght(i,:) = binavg(yday6h-0.125, TGT.yday, TGT.hght(i,:));
end

% interpolate Yap to 6 h
[Y.uwnd, Y.vwnd, Y.temp, Y.thta, Y.shum, Y.hght] = deal(NaN(np,nt));
Y.lat = Yap.lat + zeros(nt,1);
Y.lon = Yap.lon + zeros(nt,1);
Yap.mday = Yap.time+366; % correct time by 366!
for i = 1:np
    Y.uwnd(i,:) = interp1(Yap.mday, Yap.uwnd(i,:), mday6h);
    Y.vwnd(i,:) = interp1(Yap.mday, Yap.vwnd(i,:), mday6h);
    Y.temp(i,:) = interp1(Yap.mday, Yap.temp(i,:), mday6h);
    Y.thta(i,:) = interp1(Yap.mday, Yap.thta(i,:), mday6h);
    Y.shum(i,:) = interp1(Yap.mday, Yap.shum(i,:), mday6h);
    Y.hght(i,:) = interp1(Yap.mday, Yap.hght(i,:), mday6h);
end

% Palau is 6-hourly during PISTON IOP
% synchronize Palau to 6 h timestamp by nearest neighbor interp
[P.uwnd, P.vwnd, P.temp, P.thta, P.shum, P.hght] = deal(NaN(np,nt));
P.lat = Palau.lat + zeros(nt,1);
P.lon = Palau.lon + zeros(nt,1);
Palau.mday = Palau.time+366; % correct time by 366!
for i = 1:np
    P.uwnd(i,:) = interp1(Palau.mday, Palau.uwnd(i,:), mday6h, 'nearest');
    P.vwnd(i,:) = interp1(Palau.mday, Palau.vwnd(i,:), mday6h, 'nearest');
    P.temp(i,:) = interp1(Palau.mday, Palau.temp(i,:), mday6h, 'nearest');
    P.thta(i,:) = interp1(Palau.mday, Palau.thta(i,:), mday6h, 'nearest');
    P.shum(i,:) = interp1(Palau.mday, Palau.shum(i,:), mday6h, 'nearest');
    P.hght(i,:) = interp1(Palau.mday, Palau.hght(i,:), mday6h, 'nearest');
end
%{
for i = 1:np
    P.uwnd(i,:) = bindata(mday6h-0.125, Palau.mday, Palau.uwnd(i,:));
    P.vwnd(i,:) = bindata(mday6h-0.125, Palau.mday, Palau.vwnd(i,:));
    P.temp(i,:) = bindata(mday6h-0.125, Palau.mday, Palau.temp(i,:));
    P.thta(i,:) = bindata(mday6h-0.125, Palau.mday, Palau.thta(i,:));
    P.shum(i,:) = bindata(mday6h-0.125, Palau.mday, Palau.shum(i,:));
    P.hght(i,:) = bindata(mday6h-0.125, Palau.mday, Palau.hght(i,:));
end
%}

% form vectors [Palau; Yap; TGT]
% choose a reference latitude
coslat = mean(nanmean(cosd([P.lat; Y.lat; T.lat]))); % 0.985
r_earth = 6.371e6; % m

y = r_earth*([P.lat Y.lat T.lat]);        % meridional coordinate, meters
x = r_earth*coslat*([P.lon Y.lon T.lon]); % zonal coordinate, meters
u = cat(3, P.uwnd, Y.uwnd, T.uwnd);
v = cat(3, P.vwnd, Y.vwnd, T.vwnd);
% mean sounding at mid-layer
temp = cat(3, P.temp, Y.temp, T.temp); % K
hght = cat(3, P.hght, Y.hght, T.hght); % geopotential height
shum = cat(3, P.shum, Y.shum, T.shum); % kg/kg
qbar = nanmean(shum,3);
tbar = nanmean(temp,3);
ubar = nanmean(u,3);
vbar = nanmean(v,3);

Cp = 1005.7; % J/K/kg
g = 9.8;
s = Cp*temp + g*hght; % dry static energy of layers, J/kg
h = s + Lv(temp+273).*shum; % moist static energy 
sbar = nanmean(s,3);
[slev, qlev] = deal(zeros(np+1,nt));
for it=1:nt % interpolate s,q to pressure level interfaces
    qlev(:,it) = interp1(pmid,qbar(:,it), plev,'linear','extrap');
    slev(:,it) = interp1(pmid,sbar(:,it), plev,'linear','extrap');
end
su = s.*u;
sv = s.*v;
qu = shum.*u;
qv = shum.*v;
hu = h.*u;
hv = h.*v;

% centered time derivatives
dsbardt = (sbar(:,1:end-2)-sbar(:,3:end))/(2*6*3600);
dqbardt = (qbar(:,1:end-2)-qbar(:,3:end))/(2*6*3600);

% spatial derivatives from the 3 stations
[dudx, dudy, dvdx, dvdy,]                  = deal(NaN(np,nt));
[dsudx, dsvdy, dhudx, dhvdy, dqudx, dqvdy] = deal(NaN(np,nt));
[dsdx, dsdy, dqdx, dqdy]                   = deal(NaN(np,nt));
for i = 1:np
    dudx(i,:)  = partiald(x, squeeze( u(i,:,:)), 2)';
    dudy(i,:)  = partiald(y, squeeze( u(i,:,:)), 2)';
    dvdx(i,:)  = partiald(x, squeeze( v(i,:,:)), 2)';
    dvdy(i,:)  = partiald(y, squeeze( v(i,:,:)), 2)';
    % flux form derivatives
    dsudx(i,:) = partiald(x, squeeze(su(i,:,:)), 2)';
    dsvdy(i,:) = partiald(y, squeeze(sv(i,:,:)), 2)';
    dqudx(i,:) = partiald(x, squeeze(qu(i,:,:)), 2)';
    dqvdy(i,:) = partiald(y, squeeze(qv(i,:,:)), 2)';
    dhudx(i,:) = partiald(x, squeeze(hu(i,:,:)), 2)';
    dhvdy(i,:) = partiald(y, squeeze(hv(i,:,:)), 2)';
    % advective form derivatives
    dsdx(i,:)  = partiald(x, squeeze( s(i,:,:)), 2)';
    dsdy(i,:)  = partiald(x, squeeze( s(i,:,:)), 2)';
    dqdx(i,:)  = partiald(y, squeeze( shum(i,:,:)), 2)';
    dqdy(i,:)  = partiald(y, squeeze( shum(i,:,:)), 2)';
end
divg = dudx + dvdy;
vort = dvdx - dudy;
omga = [zeros(1,nt); -cumsum(100*dp.*divg)]; % Pa/s on pressure level interfaces, starting at pint(1)= 0 hPa
% does omega need some correction?
% integrating down gives large vertical velocity near surface
% adjust so that 1000 hPa omega = 0; refinements would be for Eulerian dpsfc/dt
% and convergence between the surface and p=1000hPa
omgadj = omga - bsxfun(@times, omga(end,:), plev/1000); % adjusts omga(1000hPa)=0
omgmid = (omgadj(1:end-1,:) + omgadj(2:end,:))/2;

% MSE budget terms
divsV = dsudx + dsvdy;
divqV = dqudx + dqvdy;
somg = slev.*omgadj;
qomg = qlev.*omgadj; % (np+!,nt)
dsomgdp = diff(somg,1,1)./(dp*100); 
dqomgdp = diff(qomg,1,1)./(dp*100);

Q1 =                  dsbardt(:,[1 1:end end]) + divsV + dsomgdp  ;
Q2 = -Lv(tbar+273).*( dqbardt(:,[1 1:end end]) + divqV + dqomgdp );

% Q1 and Q2 in advective form, neglecting dsbardt
ugrads = ubar.*dsdx + vbar.*dsdy;
ugradq = ubar.*dqdx + vbar.*dqdy;
dsdp = diff(slev,1,1)./dp;
dqdp = diff(qlev,1,1)./dp;
omgdsdp = omgmid.*dsdp;
omgdqdp = omgmid.*dqdp;
Q1adv =                                            ugrads + omgdsdp; % neglect dsbardt
Q2adv = -Lv(tbar+273).*(dqbardt(:,[1 1:end end]) + ugradq + omgdqdp); % USE dqbardt

% yday 251, Sept 8 ends 1st cruise

% define plot limits and ticks
xlm = [22.375 round(8*(min(datenum(2018,10,9), now-9/24)-datenum(2018,8,0)))/8];
ylm = [200 1000];
xtk = datenum(2018,8,floor(xlm(1)):ceil(xlm(2))) - datenum(2018,8,0);
[~,~,xtl] = datevec(xtk+datenum(2018,8,0));
xtlbl = num2cell(xtl);
xtlbl(mod(xtl,5)>0) = {[]};
ytk = 0:200:1000;
ytl = ytk;


% plot divergence, vorticity, and pressure velocity
b2rcolormap(17);
clf
ax(1)=subplot(3,1,1);
pcolor(mday6h([1:end end])-datenum(2018,8,0)-0.125,plev,divg([1:end end],[1:end end])); shading flat; axis ij
title('divergence (s^{-1})','fontweight','normal')
caxis([-1 1]*1e-6);
cb = colorbar; set(cb,'tickdir','out')
ylabel('pressure (hPa)')
ax(2)=subplot(3,1,2);
pcolor(mday6h([1:end end])-datenum(2018,8,0)-0.125,plev,vort([1:end end],[1:end end])); shading flat; axis ij
title('vorticity (s^{-1})','fontweight','normal')
% xlabel('August day')
caxis([-1 1]*1e-6);
cb = colorbar; set(cb,'tickdir','out')
ax(3)=subplot(3,1,3);
pcolor(mday6h([1:end end])-datenum(2018,8,0)-0.125,plev,omgadj(:,[1:end end])); shading flat; axis ij
hold on
plot(mday6h-datenum(2018,8,0),omga(end,:)*1e4+700,'k-')
plot(mday6h-datenum(2018,8,0),700+0*mday6h,'k:')
title('pressure velocity (Pa s^{-1})','fontweight','normal')
% xlabel('August day')
xlabel('          August                       September')
caxis([-1.2 1.2]*1e-2);
cb = colorbar; set(cb,'ydir','reverse','tickdir','out')
set(ax(:),'fontsize',15,'xlim',xlm,'ylim',ylm,'xtick',xtk,'xticklabel',xtlbl,'ytick',ytk,'tickdir','out')
orient tall
saveas(gcf,'divg_vort_omga_triangle.png')

% plot Q1, Q2 - looks bad from flux form
clf; ax=[];
ax(1)=subplot(3,1,1);
pcolor(mday6h([1:end end])-datenum(2018,8,0)-0.125,plev,Q1([1:end end],[1:end end])); shading flat; axis ij
title('apparent heat source Q1 (J kg^{-1} s^{-1})','fontweight','normal')
caxis([-1 1]*1e-1);
cb = colorbar; set(cb,'tickdir','out')
ylabel('pressure (hPa)')
ax(2)=subplot(3,1,2);
pcolor(mday6h([1:end end])-datenum(2018,8,0)-0.125,plev,Q2([1:end end],[1:end end])); shading flat; axis ij
title('apparent moisture sink Q2 (J kg^{-1} s^{-1})','fontweight','normal')
xlabel('August day')
caxis([-1 1]*1.e-1);
cb = colorbar; set(cb,'tickdir','out')
ax(3)=subplot(3,1,3);
pcolor(mday6h([1:end end])-datenum(2018,8,0)-0.125,plev,omgadj(:,[1:end end])); shading flat; axis ij
hold on
plot(mday6h-datenum(2018,8,0),omga(end,:)*1e4+700,'k-')
plot(mday6h-datenum(2018,8,0),700+0*mday6h,'k:')
title('pressure velocity (Pa s^{-1})','fontweight','normal')
% xlabel('August day')
xlabel('          August                       September')
caxis([-1.2 1.2]*1e-2);
cb = colorbar; set(cb,'ydir','reverse','tickdir','out')
set(ax(:),'fontsize',15,'xlim',xlm,'ylim',ylm,'xtick',xtk,...
    'xticklabel',xtlbl,'ytick',ytl,'tickdir','out')
orient tall
saveas(gcf,'Q1Q2_triangle.png')

% plot Q1 source terms. Vertical and horizontal flux divergences
% compensate.
% Q1 =                  dsbardt(:,[1 1:end end]) + divsV + dsomgdp  ;
dsbardt_ext = dsbardt(:,[1 1:end end]);
clf; ax=[];
ax(1)=subplot(3,1,1);
pcolor(mday6h([1:end end])-datenum(2018,8,0)-0.125,plev,dsbardt_ext([1:end end],[1:end end])); shading flat; axis ij
title('dsbardt (J kg^{-1} s^{-1})','fontweight','normal')
caxis([-1 1]*1e-1);
cb = colorbar; set(cb,'tickdir','out')
ylabel('pressure (hPa)')
ax(2)=subplot(3,1,2);
pcolor(mday6h([1:end end])-datenum(2018,8,0)-0.125,plev,divsV([1:end end],[1:end end])); shading flat; axis ij
title('divsV (J kg^{-1} s^{-1})','fontweight','normal')
caxis([-2 2]*1.e-1);
cb = colorbar; set(cb,'tickdir','out')
ax(3)=subplot(3,1,3);
pcolor(mday6h([1:end end])-datenum(2018,8,0)-0.125,plev,dsomgdp([1:end end],[1:end end])); shading flat; axis ij
title('dsomgdp (J kg^{-1} s^{-1})','fontweight','normal')
% xlabel('August day')
xlabel('          August                       September')
caxis([-2 2]*1.e-1);
cb = colorbar; set(cb,'tickdir','out')
set(ax(:),'fontsize',15,'xlim',xlm,'ylim',ylm,'xtick',xtk,...
    'xticklabel',xtlbl,'ytick',ytl,'tickdir','out')
orient tall
saveas(gcf,'Q1_terms.png')


% plot Q2 source terms
% Q2 = -Lv(tbar+273).*( dqbardt(:,[1 1:end end]) + divqV + dqomgdp );
Ldqbardt_ext = -Lv(tbar+273).*dqbardt(:,[1 1:end end]);
LdivqV       = -Lv(tbar+273).*divqV;
Ldqomgdp     = -Lv(tbar+273).*dqomgdp;
clf; ax=[];
ax(1)=subplot(3,1,1);
pcolor(mday6h([1:end end])-datenum(2018,8,0)-0.125,plev,Ldqbardt_ext([1:end end],[1:end end])); shading flat; axis ij
title('-Ldqbardt (J kg^{-1} s^{-1})','fontweight','normal')
caxis([-2 2]*1e-1);
cb = colorbar; set(cb,'tickdir','out')
ylabel('pressure (hPa)')
ax(2)=subplot(3,1,2);
pcolor(mday6h([1:end end])-datenum(2018,8,0)-0.125,plev,LdivqV([1:end end],[1:end end])); shading flat; axis ij
title('-LdivqV (J kg^{-1} s^{-1})','fontweight','normal')
caxis([-2 2]*1.e-2);
cb = colorbar; set(cb,'tickdir','out')
ax(3)=subplot(3,1,3);
pcolor(mday6h([1:end end])-datenum(2018,8,0)-0.125,plev,Ldqomgdp([1:end end],[1:end end])); shading flat; axis ij
title('-Ldqomgdp (J kg^{-1} s^{-1})','fontweight','normal')
% xlabel('August day')
xlabel('          August                       September')
caxis([-2 2]*1.e-2);
cb = colorbar; set(cb,'tickdir','out')
set(ax(:),'fontsize',15,'xlim',xlm,'ylim',ylm,'xtick',xtk,...
    'xticklabel',xtlbl,'ytick',ytl,'tickdir','out')
orient tall
saveas(gcf,'Q2_terms.png')

% plot Q1, Q2 in advective form; works better
clf; ax=[];
ax(1)=subplot(3,1,1);
pcolor(mday6h([1:end end])-datenum(2018,8,0)-0.125,plev,Q1adv([1:end end],[1:end end])); shading flat; axis ij
title('apparent heat source Q1 - advective form (J kg^{-1} s^{-1})','fontweight','normal')
caxis([-1 1]);
cb = colorbar; set(cb,'tickdir','out')
ylabel('pressure (hPa)')
ax(2)=subplot(3,1,2);
pcolor(mday6h([1:end end])-datenum(2018,8,0)-0.125,plev,Q2adv([1:end end],[1:end end])); shading flat; axis ij
title('apparent moisture sink Q2 - advective form (J kg^{-1} s^{-1})','fontweight','normal')
caxis([-1 1]);
cb = colorbar; set(cb,'tickdir','out')
ax(3)=subplot(3,1,3);
pcolor(mday6h([1:end end])-datenum(2018,8,0)-0.125,plev,omgadj(:,[1:end end])); shading flat; axis ij
hold on
plot(mday6h-datenum(2018,8,0),omga(end,:)*1e4+700,'k-')
plot(mday6h-datenum(2018,8,0),700+0*mday6h,'k:')
title('pressure velocity (Pa s^{-1})','fontweight','normal')
% xlabel('August day')
xlabel('          August                       September')
caxis([-1.2 1.2]*1e-2);
cb = colorbar; set(cb,'ydir','reverse','tickdir','out')
set(ax(:),'fontsize',15,'xlim',xlm,'ylim',ylm,'xtick',xtk,...
    'xticklabel',xtlbl,'ytick',ytl,'tickdir','out')
orient tall
saveas(gcf,'Q1Q2adv_triangle.png')

% plot Q1adv source terms
% Q1adv = ugrads + omgdsdp; % neglect dsbardt
dsbardt_ext = dsbardt(:,[1 1:end end]);
clf; ax=[];
ax(1)=subplot(3,1,1);
pcolor(mday6h([1:end end])-datenum(2018,8,0)-0.125,plev,dsbardt_ext([1:end end],[1:end end])); shading flat; axis ij
title('\partial{s}/\partial{t} (J kg^{-1} s^{-1})','fontweight','normal')
caxis([-1 1]*1e-1);
cb = colorbar; set(cb,'tickdir','out')
ylabel('pressure (hPa)')
ax(2)=subplot(3,1,2);
pcolor(mday6h([1:end end])-datenum(2018,8,0)-0.125,plev,ugrads([1:end end],[1:end end])); shading flat; axis ij
title('{\bf V}\cdot\nabla{s} (J kg^{-1} s^{-1})','fontweight','normal')
caxis([-2 2]*1.e-3);
cb = colorbar; set(cb,'tickdir','out')
ax(3)=subplot(3,1,3);
pcolor(mday6h([1:end end])-datenum(2018,8,0)-0.125,plev,omgdsdp([1:end end],[1:end end])); shading flat; axis ij
title('\omega \partial{s}/\partial{p} (J kg^{-1} s^{-1})','fontweight','normal')
xlabel('          August                       September')
caxis([-8 8]*1.e-1);
cb = colorbar; set(cb,'ydir','reverse','tickdir','out')
set(ax(:),'fontsize',15,'xlim',xlm,'ylim',ylm,'xtick',xtk,...
    'xticklabel',xtl,'ytick',ytl,'tickdir','out')
orient tall
saveas(gcf,'Q1adv_terms.png')

% plot Q2adv source terms
% Q2adv = -Lv(tbar+273).*(dqbardt(:,[1 1:end end]) + ugradq + omgdqdp); % USE dqbardt
Ldqbardt_ext = -Lv(tbar+273).*dqbardt(:,[1 1:end end]); % as above
Lugradq      = -Lv(tbar+273).*ugradq;
Lomgdqdp     = -Lv(tbar+273).*omgdqdp;
clf; ax=[];
ax(1)=subplot(3,1,1);
pcolor(mday6h([1:end end])-datenum(2018,8,0)-0.125,plev,Ldqbardt_ext([1:end end],[1:end end])); shading flat; axis ij
title('-L\partial{q}/\partial{t} (J kg^{-1} s^{-1})','fontweight','normal')
caxis([-2 2]*1e-1);
cb = colorbar; set(cb,'tickdir','out')
ylabel('pressure (hPa)')
ax(2)=subplot(3,1,2);
pcolor(mday6h([1:end end])-datenum(2018,8,0)-0.125,plev,Lugradq([1:end end],[1:end end])); shading flat; axis ij
title('-L{\bf{V}}\cdot\nabla{q} (J kg^{-1} s^{-1})','fontweight','normal')
caxis([-2 2]*1.e-2);
cb = colorbar; set(cb,'tickdir','out')
ax(3)=subplot(3,1,3);
pcolor(mday6h([1:end end])-datenum(2018,8,0)-0.125,plev,Lomgdqdp([1:end end],[1:end end])); shading flat; axis ij
title('-L\omega\partial{q}/\partial{p} (J kg^{-1} s^{-1})','fontweight','normal')
% xlabel('August day')
xlabel('          August                       September')
caxis([-4 4.001]*1.e-1);
cb = colorbar; set(cb,'tickdir','out')
set(ax(:),'fontsize',15,'xlim',xlm,'ylim',ylm,'xtick',xtk,...
    'xticklabel',xtl,'ytick',ytl,'tickdir','out')
orient tall
saveas(gcf,'Q2adv_terms.png')

clf
subplot(2,2,1)
plot(T.lon,T.lat,'.-')
hold on 
plot(P.lon,P.lat,'ro')
plot(Y.lon,Y.lat,'kv')
set(gca,'fontsize',12, 'dataaspectratio',[1 1 1],'xlim',[130 140])
saveas(gcf,'sounding_array_plan.eps')

% quick plot
pcolor(mday6h, pmid, T.uwnd); shading flat; datetick; axis ij; caxis([-30 30]); axis tight
pcolor(mday6h, pmid, h); shading flat; datetick; axis ij; caxis([-30 30]); axis tight