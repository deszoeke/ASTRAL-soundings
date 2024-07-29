cd('/Users/sdeszoek/Data/cruises/PISTON_MISOBOB_2018/Soundings/Thompson/matlab');

addpath(genpath('/Users/sdeszoek/Documents/MATLAB/user_tools/stat/'));
addpath(genpath('/Users/sdeszoek/Documents/MATLAB/user_tools/thermo/'));
addpath(genpath('/Users/sdeszoek/Documents/MATLAB/user_tools/graphics/'));
addpath(genpath('/Users/sdeszoek/Documents/MATLAB/user_tools/fileio/'));

filepath = '/Users/sdeszoek/Data/cruises/PISTON_MISOBOB_2018/Soundings/Thompson/netcdf';
flist = dir(fullfile(filepath,'Thompson*.nc'));
nsonde = length(flist);


%timestamp = sprintf('%4d%02d%02d_%02d%02d%02d', year, month, day, hour, minute, second);
%ncf = fullfile(filepath, ['Thompson_sonde' timestamp '.nc']);

zint = (0:10:27000)';
nz = length(zint);
[pres, temp, rhum, uwnd, vwnd, shum, thta, thte, dwpt] = deal(NaN(nz, nsonde));
[time, lat, lon] = deal(NaN(1,nsonde));
for i = 1:nsonde
    ncf = fullfile(filepath, flist(i).name);
    
    sonde = read_TGT_sonde(ncf);
    time(i) = nanmean(sonde.Time); % seconds since Jan 1, 2018 = (yday-1)*86400
    lat(i) = nanmean(sonde.Latitude_of_balloon); % seconds since Jan 1, 2018 = (yday-1)*86400
    lon(i) = nanmean(sonde.Longitude_of_balloon); % seconds since Jan 1, 2018 = (yday-1)*86400
    [t1, t2, t3, t4, t5, t6, t7, t8, t9] = interp_sonde_z(sonde, zint);
    pres(:,i) = t1;
    temp(:,i) = t2; % K
    rhum(:,i) = t3;
    uwnd(:,i) = t4;
    vwnd(:,i) = t5;
    shum(:,i) = t6;
    thta(:,i) = t7;
    thte(:,i) = t8;
    dwpt(:,i) = t9;
end

% Pi = temp./thta;
thdw = dwpt.*(thta./temp); % dew point potential temperature
thes = theta_es(pres*100, temp); % saturation equivalent potential temperature
qast = 1.0e3*qs(pres*100,temp-273.15); % saturation specific humidity

thtw = theta_w(pres*100, temp, shum/1e3); % interesting how constant over troposphere
wspd = uwnd.^2 + vwnd.^2;

yday = time/86400 + 1;
% For some reason 20180906_173022 (249.7687) is a duplicate of 20180906_143022.
% A look at the edt files shows no 17:30:22 sounding.
% 20180906_173022 is removed by renaming it to 
% duplicates_20180906_143022_Thompson_sonde_20180906_173022.nc
% and the launch time for yday(126) should be 14:30:22 (249.6044).
if yday(126) > 249.75 % hack until fixed in netcdf files
    yday(126) = datenum(2018,09,06,14,30,22)-datenum(2018,0,0);
end

anom = @(x,d) x - nanmean(x,d);
plot(anom(thta,2),zint)

clf
ax=zeros(4,1);
ax(1)=subplot(2,2,1);
plot(nanmean(uwnd,2),zint/1e3,'k-','linewidth',3); title('u (m/s)')
hold on
plot(uwnd,zint/1e3,'linewidth',0.5)
ax(2)=subplot(2,2,2);
plot(nanmean(vwnd,2),zint/1e3,'k-','linewidth',3); title('v (m/s)')
hold on
plot(vwnd,zint/1e3,'linewidth',0.5)

ax(3)=subplot(2,2,3);
plot(thdw,zint/1e3,'b','linewidth',0.5); hold on
plot(thte,zint/1e3,'r','linewidth',0.5)
plot(thta,zint/1e3,'k','linewidth',0.5)
plot(thtw,zint/1e3,'color',[0 0.7 0.1],'linewidth',0.5)
% plot(nanmean(thdw,2),zint/1e3,'b-','linewidth',3);
% plot(nanmean(thte,2),zint/1e3,'r-','linewidth',3);
% plot(nanmean(thta,2),zint/1e3,'k-','linewidth',3); title('theta (K)')
xlim([290, 360])
ax(4)=subplot(2,2,4);
plot(shum,zint/1e3,'linewidth',0.5); hold on
plot(qast,zint/1e3,'r','linewidth',0.5);
% plot(nanmean(shum,2),zint/1e3,'-','linewidth',3,'color',0.7*[1 1 1]); title('qv (g/kg)')
% plot(nanmean(qast,2),zint/1e3,'r-','linewidth',3);
set(ax(4),'xscale','log','xlim',[1e-3 50])
set(ax(:),'fontsize',14,'ylim',[0 18])
% saveas(gcf,'TGT_sondes_4pane.eps','epsc')

clf
plot(thdw,zint/1e3,'b','linewidth',0.5); hold on
plot(thte,zint/1e3,'r','linewidth',0.5)
plot(thta,zint/1e3,'k','linewidth',0.5)
plot(thtw,zint/1e3,'color',[0 0.7 0.1],'linewidth',0.5)
% plot(nanmean(thdw,2),zint/1e3,'b-','linewidth',3);
% plot(nanmean(thte,2),zint/1e3,'r-','linewidth',3);
% plot(nanmean(thta,2),zint/1e3,'k-','linewidth',3); title('theta (K)')
xlim([290, 360]); ylim([0 16])
set(gca,'fontsize',16)
% legend('dewpoint','equivalent','\theta','wet bulb'); legend boxoff
text(300,15,'wet bulb','fontsize',16)
text(322,15,'dew point','fontsize',16)
text(314,2 ,'\theta','fontsize',16)
text(348,5 ,'equivalent','fontsize',16)
xlabel('potential temperature (K)')
ylabel('height (km)')
% print('-dpng','thetas.png')

% define plot limits and ticks
xlm = [234 round(8*(min(datenum(2018,0,281), now)-datenum(2018,0,0)))/8];
ylm = [0 20];
xtk = datenum(2018,0,floor(xlm(1)):1:ceil(xlm(2))) - datenum(2018,0,0);
[~,~,xtl] = datevec(xtk+datenum(2018,0,0));
xtlbl = num2cell(xtl);
xtlbl(mod(xtl,5)>0) = {[]};
ytk = 0:5:20;
ytl = num2cell(ytk);
ytl(mod(ytk,5)>0) = {[]};

clf; ax=[];
b2rcolormap(33);
ax(1)=subplot(2,1,1);
pcolor(yday,zint/1e3,thta); shading flat
colorbar
ax(2)=subplot(2,1,2);
pcolor(yday,zint/1e3,shum); shading flat
colorbar
set(ax(:),'fontsize',16,'xlim',xlm,'ylim',ylm,'ytick',ytk,'yticklabel',ytl,'xtick',xtk,'xticklabel',xtlbl,'ytick',ytk,'tickdir','out')

clf; ax=[];
b2rcolormap(25);
ax(1)=subplot(2,1,1);
pcolor(yday,zint/1e3,anom(uwnd,2)); shading flat
caxis([-12 12])
colorbar
ylabel('height (km)')
title('zonal wind anomaly (m/s)','fontweight','normal')
ax(2)=subplot(2,1,2);
pcolor(yday,zint/1e3,anom(vwnd,2)); shading flat
caxis([-12 12])
xlabel('August                               September                                     October')
title('meridional wind anomaly (m/s)','fontweight','normal')
colorbar
set(ax(:),'fontsize',16,'xlim',xlm,'ylim',ylm,'ytick',ytk,'yticklabel',ytl,'xtick',xtk,'xticklabel',xtlbl,'ytick',ytk,'tickdir','out')
set(ax(:),'yminortick','on')
orient landscape
print('-dpng','sonde_anom_uv_timeheight.png','-r300')

clf; ax=[];
b2rcolormap(17);
ax(1)=subplot(2,1,1);
pcolor(yday,zint/1e3,anom(thta,2)); shading flat
caxis([-3 3])
colorbar
ylabel('height (km)')
title('potential temperature anomaly (K)','fontweight','normal')
ax(2)=subplot(2,1,2);
pcolor(yday,zint/1e3,anom(shum,2)); shading flat
title('specific humidity anomaly (g/kg)','fontweight','normal')
caxis([-2 2])
% pcolor(yday,zint/1e3,anom(shum,2)./nanmean(shum,2)); shading flat
% title('specific humidity anomaly (fraction of mean)','fontweight','normal')
% caxis([-1 1])
xlabel('yearday')
colorbar
set(ax(:),'fontsize',16,'xlim',xlm,'ylim',ylm,'xtick',xtk,'xticklabel',xtlbl,'ytick',ytk,'tickdir','out')
set(ax(:),'ylim',[0 12])
set(gca,'yminortick','on')
print('-dpng','sonde_anom_th_q_timeheight.png','-r300')

% scatter plot
% thtw and q at 2 km are 
clf
plot(thtw(201,:),shum(201,:),'.') % 2 km
xlabel('\theta_w (K)'); ylabel('q (g/kg)')
hold on

% diurnal cycle
% diurnal cycle below 2 km
dielbin = (-1.5:3:25.5)'/24; % cyclical bins
dielbinctr = (0:3:24)';
nbin = length(dielbinctr)-1;
dieli = interp1(dielbinctr,1:length(dielbinctr),mod(24*yday,24),'nearest');
dieli(dieli==9) = 1; % wrap last bin 24Z to 00Z

% accumulate
[dithta,diwspd,dishum] = deal(zeros(nz, nbin));
[nthta,nwspd,nshum]    = deal(zeros(nz, nbin));
for i = 1:nsonde
    kk = isfinite(thta(:,i));
    dithta(kk,dieli(i)) = dithta(kk,dieli(i)) + thta(kk,i);
    nthta( kk,dieli(i)) = nthta(kk,dieli(i)) + 1;
    kk = isfinite(wspd(:,i));
    diwspd(kk,dieli(i)) = diwspd(kk,dieli(i)) + wspd(kk,i);
    nwspd( kk,dieli(i)) = nwspd(kk,dieli(i)) + 1;
    kk = isfinite(shum(:,i));   
    dishum(kk,dieli(i)) = dishum(kk,dieli(i)) + shum(kk,i);
    nshum( kk,dieli(i)) = nshum(kk,dieli(i)) + 1;
end
dithta = dithta./nthta;
diwspd = diwspd./nwspd;
dishum = dishum./nshum;

clf
b2rcolormap(9);
contourf(0:3:30,zint(1:331)/1e3,anom(dithta(1:331,[1:end 1:3]),2),-0.4:0.1:0.4,'edgecolor','none');
colorbar
set(gca,'xtick',0:3:30, 'xticklabel',[0:3:24 0:3:6], 'fontsize',16)
set(gca,'xticklabel',[0:3:21 0:3:6],'fontsize',16)
ylabel('height (m)'); xlabel('UTC hour')
hold on
contour(0:3:30, zint(1:331)/1e3, dithta(1:331,[1:end 1:3]), 280:320,'k');
contour(0:3:30, zint(1:331)/1e3, anom(dishum(1:331,[1:end 1:3]),2), -0.5:0.1:-0.1, 'color',[0.7 .2 0],'linewidth',1.5); % reds
contour(0:3:30, zint(1:331)/1e3, anom(dishum(1:331,[1:end 1:3]),2), 0.1:0.1:0.5, 'color',[0.3 0.7 .3],'linewidth',1.5); % greens
caxis([-0.4, 0.4])
plot(0:0.25:30, max(0,cos(2*pi*((0:0.25:30)+9+12)/24)),'k-','linewidth',2)
saveas(gcf,'diurnal_0-3km_all.png')

% cold days might have a different diurnal cycle
datenum(2018,9,[2 5])-datenum(2018,0,0);
icold = datenum(2018,9,2)-datenum(2018,0,0) <= floor(yday) <= datenum(2018,9,5)-datenum(2018,0,0);

%{
% Try computing CAPE in a weird way using theta_e_virt of the parcel
dz = 10; % m
g = 9.8; % m/s
basez = 700;
virt = 1.0 + 0.61*shum/1e3;
base = find(zint>=basez,1,'first');
thev_parcel = nanmean(thte(1:base,:).*virt(1:base,:)); % mix out between 0-700 m
thtv = thta.*virt;
% % try comparing theta_e_v_parcel with theta_es_v_environment
top = find(bsxfun(@and, zint>base & zint<20e3, thev_parcel<thtv), 1, 'first');
cape = dz*g*sum( max(0, thev_parcel./thtv(base:top) - 1.0) );
%                       !!!!!!!!!!!
% This CAPE calc is wrong. Either integrate T_parcel up from Tlcl, or invert the
% transcendental equation for theta_e of the moist adiabat.
%}
