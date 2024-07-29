cd('/Users/sdeszoek/Data/cruises/PISTON_MISOBOB_2018/Soundings/Thompson/matlab')
addpath(genpath('/Users/sdeszoek/Documents/MATLAB/user_tools/stat/'));
addpath(genpath('/Users/sdeszoek/Documents/MATLAB/user_tools/thermo/'));
addpath(genpath('/Users/sdeszoek/Documents/MATLAB/user_tools/graphics/'));
addpath(genpath('/Users/sdeszoek/Documents/MATLAB/user_tools/fileio/'));

clear dir z
filepath = '/Users/sdeszoek/Data/cruises/PISTON_MISOBOB_2018/Soundings/Thompson/netcdf';
flist = dir(fullfile(filepath,'Thompson*.nc'));
nsonde = length(flist);

%timestamp = sprintf('%4d%02d%02d_%02d%02d%02d', year, month, day, hour, minute, second);
%ncf = fullfile(filepath, ['Thompson_sonde' timestamp '.nc']);

ncrwb = fullfile('/Users/sdeszoek/Data/cruises/PISTON_MISOBOB_2018/TGT/Wband/Mom_4/20182450000MMCRMom.nc');
zint = 11.7 + ncread(ncrwb,'Heights',[1 1],[Inf 1]); % m, start (bottom) of gates
% first valid gate is zint(3) = 108 m
nz = length(zint);
zint(nz+1) = zint(nz) + zint(nz)-zint(nz-1);
[pres, temp, rhum, uwnd, vwnd, shum, thta, thte, dwpt, wspd, wdir] = deal(NaN(nz, nsonde));
[time, lat, lon] = deal(NaN(1,nsonde));
for i = 1:nsonde
    ncf = fullfile(filepath, flist(i).name);
    
    sonde(i) = read_TGT_sonde(ncf);
    time(i) = nanmean(sonde(i).Time); % seconds since Jan 1, 2018 = (yday-1)*86400
    lat(i) = nanmean(sonde(i).Latitude_of_balloon); % seconds since Jan 1, 2018 = (yday-1)*86400
    lon(i) = nanmean(sonde(i).Longitude_of_balloon); % seconds since Jan 1, 2018 = (yday-1)*86400
    [t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11] = avg_sonde_z_hed(sonde(i), zint);
    pres(:,i) = t1(1:end-1);
    temp(:,i) = t2(1:end-1); % K
    rhum(:,i) = t3(1:end-1);
    uwnd(:,i) = t4(1:end-1);
    vwnd(:,i) = t5(1:end-1);
    shum(:,i) = t6(1:end-1);
    thta(:,i) = t7(1:end-1);
    thte(:,i) = t8(1:end-1);
    dwpt(:,i) = t9(1:end-1);
    wspd(:,i) = t10(1:end-1);
    wdir(:,i) = t11(1:end-1);
end

% Pi = temp./thta;
thdw = dwpt.*(thta./temp); % dew point potential temperature
thes = theta_es(pres*100, temp); % saturation equivalent potential temperature
qast = 1.0e3*qs(pres*100,temp-273.15); % saturation specific humidity

thtw = theta_w(pres*100, temp, shum/1e3); % interesting how constant over troposphere

yday = time/86400 + 1;

% For some reason 20180906_173022 (249.7687) is a duplicate of 20180906_143022.
% A look at the edt files shows no 17:30:22 sounding.
% 20180906_173022 is removed by renaming it to 
% duplicates_20180906_143022_Thompson_sonde_20180906_173022.nc
% and the launch time for yday(126) should be 14:30:22 (249.6044).
if yday(126) > 249.75 % hack until fixed in netcdf files
    yday(126) = datenum(2018,09,06,14,30,22)-datenum(2018,0,0);
end

save('sondes_4_wband.mat',...
     'time','lat','lon','yday','zint',...
     'pres','temp','rhum','uwnd','vwnd',...
     'shum','thta','thte','dwpt','wspd',...
     'wdir','thdw','thes','qast','thtw'    );

zmid = (zint(1:end-1)+zint(2:end))/2;
anom = @(x,d) bsxfun(@minus, x, nanmean(x,d));
plot(anom(thta,2),zmid)

clf
ax=zeros(4,1);
ax(1)=subplot(2,2,1);
plot(nanmean(uwnd,2),zmid/1e3,'k-','linewidth',3); title('u (m/s)')
hold on
plot(uwnd,zmid/1e3,'linewidth',0.5)
ax(2)=subplot(2,2,2);
plot(nanmean(vwnd,2),zmid/1e3,'k-','linewidth',3); title('v (m/s)')
hold on
plot(vwnd,zmid/1e3,'linewidth',0.5)

ax(3)=subplot(2,2,3);
plot(thdw,zmid/1e3,'b','linewidth',0.5); hold on
plot(thte,zmid/1e3,'r','linewidth',0.5)
plot(thta,zmid/1e3,'k','linewidth',0.5)
plot(thtw,zmid/1e3,'color',[0 0.7 0.1],'linewidth',0.5)
% plot(nanmean(thdw,2),zmid/1e3,'b-','linewidth',3);
% plot(nanmean(thte,2),zmid/1e3,'r-','linewidth',3);
% plot(nanmean(thta,2),zmid/1e3,'k-','linewidth',3); title('theta (K)')
xlim([290, 360])
ax(4)=subplot(2,2,4);
plot(shum,zmid/1e3,'linewidth',0.5); hold on
plot(qast,zmid/1e3,'r','linewidth',0.5);
% plot(nanmean(shum,2),zmid/1e3,'-','linewidth',3,'color',0.7*[1 1 1]); title('qv (g/kg)')
% plot(nanmean(qast,2),zmid/1e3,'r-','linewidth',3);
set(ax(4),'xscale','log','xlim',[1e-3 50])
set(ax(:),'fontsize',14,'ylim',[0 18])
saveas(gcf,'TGT_sondes_4pane.eps','epsc')

clf
b2rcolormap(33)
subplot(2,1,1);
pcolor(yday,zint,thta([1:end end],:)); shading flat
colorbar
subplot(2,1,2);
pcolor(yday,zint,shum([1:end end],:)); shading flat
colorbar

clf; ax=[];
b2rcolormap(25);
ax(1)=subplot(2,1,1);
pcolor(yday,zint/1e3,anom(uwnd([1:end end],:),2)); shading flat
caxis([-12 12])
colorbar
ylabel('height (km)')
title('zonal wind anomaly (m/s)','fontweight','normal')
ax(2)=subplot(2,1,2);
pcolor(yday,zint/1e3,anom(vwnd([1:end end],:),2)); shading flat
caxis([-12 12])
xlabel('yearday')
title('meridional wind anomaly (m/s)','fontweight','normal')
colorbar
set(ax(:),'fontsize',16,'xtick',232:240,'tickdir','out')
print('-dpng','sonde_anom_uv_timeheight.png','-r300')

clf; ax=[];
b2rcolormap(25);
ax(1)=subplot(2,1,1);
pcolor(yday,zint(1:end-1)/1e3,anom(thta,2)); shading flat
caxis([-4 4])
colorbar
ylabel('height (km)')
title('potential temperature anomaly (K)','fontweight','normal')
ax(2)=subplot(2,1,2);
pcolor(yday,zint(1:end-1)/1e3,anom(shum,2)); shading flat
title('specific humidity anomaly (g/kg)','fontweight','normal')
caxis([-4 4])
% pcolor(yday,zint(1:end-1)/1e3,anom(shum,2)./nanmean(shum,2)); shading flat
% title('specific humidity anomaly (fraction of mean)','fontweight','normal')
% caxis([-1 1])
xlabel('yearday')
colorbar
xlm = [floor(min(yday)) ceil(max(yday))];
set(ax(:),'fontsize',16,'xtick',xlm(1):2:xlm(2),'xlim',xlm,'tickdir','out')

% scatter plot
clf
plot(thtw(201,:),shum(201,:),'.') % 2 km
hold on

% plot density sorting
plot(sonde(1).Potential_Temp,sonde(1).Height,'.')
hold on
plot(200 + 100*log10(sonde(1).dth.^2)/2, sonde(1).Height)

clf
plot(100*log10(sonde(1).dthe.^2)/2, sonde(1).Height)
hold on
plot(sonde(1).Equiv_Pot_Temp,sonde(1).Height,'.')
plot(sonde(1).thesort,sonde(1).Height,'.')
xlim([0 500]); ylim([0 18e3])
% one type of CAPE is ConvPE - ConvPE_groundstate 
% = integrate( g*z*(rho(theta_ev) - rho(thevsort)) dz )
% rho g dz = -dp; p/(Rd*Tv) = rho
%

%{
WRONG
% Try computing CAPE in a weird way using theta_e_virt of the parcel
dz = 10; % m
g = 9.8; % m/s
basez = 700;
virt = 1.0 + 0.61*shum/1e3;
base = find(zint>=basez,1,'first');
thev_parcel = nanmean(thte(1:base,:).*virt(1:base,:)); % mix out between 0-700 m
thtv = thta.*virt;
top = find(bsxfun(@and, zint>base & zint<20e3, thev_parcel<thtv), 1, 'first');
cape = dz*g*sum( max(0, thev_parcel./thtv(base:top) - 1.0) );
%                       !!!!!!!!!!!
% This CAPE calc is wrong. Either integrate T_parcel up from Tlcl, or invert the
% transcendental equation for theta_e of the moist adiabat.
% Maybe compare theta_e_virt of parcel and with theta_e_Sat_virt of environment???
%}
