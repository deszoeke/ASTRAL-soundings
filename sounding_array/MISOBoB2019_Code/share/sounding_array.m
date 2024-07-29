cd('/Volumes/cruise/SR1911/share/scripts/radiosonde/')

addpath(genpath('/Users/sdeszoek/Documents/MATLAB/user_tools/stat/'));
addpath(genpath('/Users/sdeszoek/Documents/MATLAB/user_tools/graphics/'));
addpath(genpath('/Users/sdeszoek/Documents/MATLAB/user_tools/thermo/'));
thermo_constants;

paren = @(X, varargin) X(varargin{:});

update = false;
if update
    % gets U.Wyo soundings:
    fetch_uwyo = 'bash /Volumes/cruise/SR1911/share/data/radiosonde/download_soundings.bsh';
    % 'curl -s "http://weather.uwyo.edu/cgi-bin/sounding?region=pac&TYPE=TEXT%3ALIST&YEAR=2018&MONTH=09&FROM=0100&TO=3012&STNM=91408" | sed -e "s/<[a-zA-Z\/][^>]*>//g" > ~/Data/cruises/PISTON_MISOBOB_2018/Soundings/Palau/91408_KororPalau/201809.txt &';
    % copy data from cruise share to laptop
    rsync_ride = 'rsync -parv /Volumes/cruise/SR1911/share/data/radiosonde/ ~/Data/cruises/MISOBOB_2019/radiosonde/';
    
    system(fetch_uwyo)
    system(rsync_ride)
    
    % read and grid data into .mat files
    ReadSondes
    ReadWyoSondes
    GridSondes
end

plev = [0:25:150 200:50:1000]';
dp = diff(plev);
pmid = plev(1:end-1) + dp/2;
np = length(pmid);

% load gridded soundings
datadir = '/Volumes/cruise/SR1911/share/data/radiosonde/gridded/';
stations = {'Bhubaneswar', 'Chennai', 'Kolkata', 'Port Blair', 'Sally Ride', 'Thiruvananthapuram', 'Vishakhapatnam' };
wordsin = @(x) x(regexp(x, '\w')); % removes spaces
Stn = cellfun(wordsin, stations, 'UniformOutput', false); % spaces removed
latlon = [ 20.2961, 85.8245; % Bhubaneswar
           13.0827, 80.2707; % Chennai
           22.5726, 88.3639; % Kolkata
           11.6234, 92.7265; % Port Blair
              15.5, 89;      % Sally Ride, but ship lat and lon vary
            8.5241, 76.9366; % Thiruvananthapuram
           17.6868, 83.2185; ]; % Vishakhapatnam

for ip = 1:length(Stn)
    A(ip) = getfield( load(fullfile(datadir, [Stn{ip} '_grd.mat'])), 'grd');
end
np0 = paren(size(A(1).temp, 1));

% 6 h time grid;
year = 2019;
mday6h = datenum(year,7,10):0.25:min(now,datenum(year,8,5));
yday6h = mday6h - datenum(year,0,0);
nt = length(yday6h);

[B.uwnd, B.vwnd, B.temp, B.thta, B.shum, B.hght] = deal(NaN(np,nt));
%interp all stations to 6h in time
for ist = 1:length(Stn)
    npgrd = size(A(ist).u,1);
    B(ist).uwnd = interp2((1:npgrd)', A(ist).time, A(ist).u',    (1:npgrd)', mday6h)';
    B(ist).vwnd = interp2((1:npgrd)', A(ist).time, A(ist).v',    (1:npgrd)', mday6h)';
    B(ist).temp = interp2((1:npgrd)', A(ist).time, A(ist).temp', (1:npgrd)', mday6h)';
    B(ist).thta = interp2((1:npgrd)', A(ist).time, A(ist).ptmp', (1:npgrd)', mday6h)';
    B(ist).shum = interp2((1:npgrd)', A(ist).time, A(ist).qv',   (1:npgrd)', mday6h)';
    gp = A(ist).gp; % use either height or geopotential
    if all(isnan(gp))
        B(ist).hght = interp2((1:npgrd)', A(ist).time, A(ist).h', (1:npgrd)', mday6h)';
    else
        B(ist).hght = interp2((1:npgrd)', A(ist).time, gp', (1:npgrd)', mday6h)';
    end
end

% choose a reference latitude
lat = latlon(:,1);
lon = latlon(:,2);
coslat = cosd(lat); % ~0.96
sinlat = sind(lat);
r_earth = 6.371e6; % m

y = r_earth * lat;          % meridional coordinate, meters
x = r_earth * coslat .* lon; % zonal coordinate, meters
% mean sounding at mid-layer
% clear u v temp hght
for istn = 1:length(Stn)
    for it = 1:nt
        u(:,it,istn)    = binavg( plev, A(1).p, B(istn).uwnd(:,it) ); % time, height
        v(:,it,istn)    = binavg( plev, A(1).p, B(istn).vwnd(:,it) );
        temp(:,it,istn) = binavg( plev, A(1).p, B(istn).temp(:,it) ); % K
        hght(:,it,istn) = binavg( plev, A(1).p, B(istn).hght(:,it) ); % geopotential height
        shum(:,it,istn) = binavg( plev, A(1).p, B(istn).shum(:,it)*1e-3 ); % kg/kg
    end
end
u = u(1:end-1,:,:);
v = v(1:end-1,:,:);
temp = temp(1:end-1,:,:);
hght = hght(1:end-1,:,:);
shum = shum(1:end-1,:,:);
qbar = nanmean(shum,3);
tbar = nanmean(temp,3);
ubar = nanmean(u,3);
vbar = nanmean(v,3);

% compute dry and moist static energies
Cp = 1005.7; % J/K/kg
g = 9.8;
s = Cp*temp + g*hght; % dry static energy of layers, J/kg
h = s + Lv(temp+273).*shum; % moist static energy 
sbar = nanmean(s,3);
[slev, qlev] = deal(zeros(np+1,nt));
for it=1:nt % interpolate s,q to pressure level interfaces
    qlev(:,it) = interp1(pmid(:),qbar(:,it), plev,'linear','extrap');
    slev(:,it) = interp1(pmid(:),sbar(:,it), plev,'linear','extrap');
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
for ip = 1:np
    dudx(ip,:)  = partiald(x', squeeze( u(ip,:,:)), 2)';
    dudy(ip,:)  = partiald(y', squeeze( u(ip,:,:)), 2)';
    dvdx(ip,:)  = partiald(x', squeeze( v(ip,:,:)), 2)';
    dvdy(ip,:)  = partiald(y', squeeze( v(ip,:,:)), 2)';
    % flux form derivatives
    dsudx(ip,:) = partiald(x', squeeze(su(ip,:,:)), 2)';
    dsvdy(ip,:) = partiald(y', squeeze(sv(ip,:,:)), 2)';
    dqudx(ip,:) = partiald(x', squeeze(qu(ip,:,:)), 2)';
    dqvdy(ip,:) = partiald(y', squeeze(qv(ip,:,:)), 2)';
    dhudx(ip,:) = partiald(x', squeeze(hu(ip,:,:)), 2)';
    dhvdy(ip,:) = partiald(y', squeeze(hv(ip,:,:)), 2)';
    % advective form derivatives
    dsdx(ip,:)  = partiald(x', squeeze( s(ip,:,:)), 2)';
    dsdy(ip,:)  = partiald(x', squeeze( s(ip,:,:)), 2)';
    dqdx(ip,:)  = partiald(y', squeeze( shum(ip,:,:)), 2)';
    dqdy(ip,:)  = partiald(y', squeeze( shum(ip,:,:)), 2)';
end
divg = dudx + dvdy;
vort = dvdx - dudy;

% adjust mean column divergence to be 0 to ensure 
% column-integrated nondivergence and omega(sfc) = omega(p=0) = 0
divg_integrand = divg - nanmean(divg(5:end,:), 1);
% Integrate omega up from zero at the surface
% Pa/s on pressure level interfaces, starting at pint(1)= 0 hPa
omga = flipud([zeros(1,nt); cumsum(flipud( 100*dp .* divg_integrand ))]);

% Don't do this; achieved zero omega at top, sfc by adjusting divergence
% omgadj = omga - bsxfun(@times, omga(end,:), plev/1000); % adjusts omga(1000hPa)=0
omgmid = (omga(1:end-1,:) + omga(2:end,:))/2;

% MSE budget terms
divsV = dsudx + dsvdy;
divqV = dqudx + dqvdy;
somg = slev.*omga;
qomg = qlev.*omga; % (np+!,nt)
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
date0 = datenum(2019,7,0);
xlm = [datenum(2019,7,10) min(ceil(now), datenum(2019,8,5))] - date0; % date of July 2019
ylm = [200 1000];
xtk = datenum(2019,7,floor(xlm(1)):ceil(xlm(2))) - date0;
[~,~,xtl] = datevec(xtk+date0);
xtlbl = num2cell(xtl);
xtlbl(mod(xtl,2)>0) = {[]};
ytk = 0:200:1000;
ytl = ytk;

% plot divergence, vorticity, and pressure velocity
b2rcolormap(17);
clf
ax(1)=subplot(3,1,1);
pcolor(mday6h([1:end end])-date0-0.25,plev,divg([1:end end],[1:end end])); shading flat; axis ij
title('divergence (s^{-1})','fontweight','normal')
caxis([-.5 .5]*1e-6);
cb = colorbar; set(cb,'tickdir','out')
ylabel('pressure (hPa)')
ax(2)=subplot(3,1,2);
pcolor(mday6h([1:end end])-date0-0.125,plev,vort([1:end end],[1:end end])); shading flat; axis ij
title('vorticity (s^{-1})','fontweight','normal')
% xlabel('August day')
caxis([-1 1]*1e-6);
cb = colorbar; set(cb,'tickdir','out')
ax(3)=subplot(3,1,3);
pcolor(mday6h([1:end end])-date0-0.25,plev,omga(:,[1:end end])); shading flat; axis ij
hold on
plot(mday6h-date0,omga(end,:)*1e4+700,'k-')
plot(mday6h-date0,700+0*mday6h,'k:')
title('pressure velocity (Pa s^{-1})','fontweight','normal')
% xlabel('August day')
xlabel('          August                       September')
caxis([-1.2 1.2]*1e-2);
cb = colorbar; set(cb,'ydir','reverse','tickdir','out')
set(ax(:),'fontsize',15,'xlim',xlm,'ylim',ylm,'xtick',xtk,'xticklabel',xtlbl,'ytick',ytk,'tickdir','out')
set(ax(:),'xlim',[10 22.5])
orient tall
saveas(gcf,'divg_vort_omga_triangle.png')

% plot Q1, Q2 - looks bad from flux form
clf; ax=[];
ax(1)=subplot(3,1,1);
pcolor(mday6h([1:end end])-date0-0.25,plev,Q1([1:end end],[1:end end])); shading flat; axis ij
title('apparent heat source Q1 (J kg^{-1} s^{-1})','fontweight','normal')
caxis([-1 1]*1e-1);
cb = colorbar; set(cb,'tickdir','out')
ylabel('pressure (hPa)')
ax(2)=subplot(3,1,2);
pcolor(mday6h([1:end end])-date0-0.25,plev,Q2([1:end end],[1:end end])); shading flat; axis ij
title('apparent moisture sink Q2 (J kg^{-1} s^{-1})','fontweight','normal')
xlabel('August day')
caxis([-1 1]*1.e-1);
cb = colorbar; set(cb,'tickdir','out')
ax(3)=subplot(3,1,3);
pcolor(mday6h([1:end end])-date0-0.25,plev,omga(:,[1:end end])); shading flat; axis ij
hold on
plot(mday6h-date0,omga(end,:)*1e4+700,'k-')
plot(mday6h-date0,700+0*mday6h,'k:')
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
pcolor(mday6h([1:end end])-date0-0.25,plev,dsbardt_ext([1:end end],[1:end end])); shading flat; axis ij
title('dsbardt (J kg^{-1} s^{-1})','fontweight','normal')
caxis([-1 1]*1e-1);
cb = colorbar; set(cb,'tickdir','out')
ylabel('pressure (hPa)')
ax(2)=subplot(3,1,2);
pcolor(mday6h([1:end end])-date0-0.25,plev,divsV([1:end end],[1:end end])); shading flat; axis ij
title('divsV (J kg^{-1} s^{-1})','fontweight','normal')
caxis([-2 2]*1.e-1);
cb = colorbar; set(cb,'tickdir','out')
ax(3)=subplot(3,1,3);
pcolor(mday6h([1:end end])-date0-0.25,plev,dsomgdp([1:end end],[1:end end])); shading flat; axis ij
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
pcolor(mday6h([1:end end])-date0-0.25,plev,Ldqbardt_ext([1:end end],[1:end end])); shading flat; axis ij
title('-Ldqbardt (J kg^{-1} s^{-1})','fontweight','normal')
caxis([-2 2]*1e-1);
cb = colorbar; set(cb,'tickdir','out')
ylabel('pressure (hPa)')
ax(2)=subplot(3,1,2);
pcolor(mday6h([1:end end])-date0-0.25,plev,LdivqV([1:end end],[1:end end])); shading flat; axis ij
title('-LdivqV (J kg^{-1} s^{-1})','fontweight','normal')
caxis([-2 2]*1.e-2);
cb = colorbar; set(cb,'tickdir','out')
ax(3)=subplot(3,1,3);
pcolor(mday6h([1:end end])-date0-0.25,plev,Ldqomgdp([1:end end],[1:end end])); shading flat; axis ij
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
pcolor(mday6h([1:end end])-date0-0.25,plev,Q1adv([1:end end],[1:end end])); shading flat; axis ij
title('apparent heat source Q1 - advective form (J kg^{-1} s^{-1})','fontweight','normal')
caxis([-1 1]);
cb = colorbar; set(cb,'tickdir','out')
ylabel('pressure (hPa)')
ax(2)=subplot(3,1,2);
pcolor(mday6h([1:end end])-date0-0.25,plev,Q2adv([1:end end],[1:end end])); shading flat; axis ij
title('apparent moisture sink Q2 - advective form (J kg^{-1} s^{-1})','fontweight','normal')
caxis([-1 1]);
cb = colorbar; set(cb,'tickdir','out')
ax(3)=subplot(3,1,3);
pcolor(mday6h([1:end end])-date0-0.25,plev,omga(:,[1:end end])); shading flat; axis ij
hold on
plot(mday6h-date0,omga(end,:)*1e4+700,'k-')
plot(mday6h-date0,700+0*mday6h,'k:')
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
pcolor(mday6h([1:end end])-date0-0.25,plev,dsbardt_ext([1:end end],[1:end end])); shading flat; axis ij
title('\partial{s}/\partial{t} (J kg^{-1} s^{-1})','fontweight','normal')
caxis([-1 1]*1e-1);
cb = colorbar; set(cb,'tickdir','out')
ylabel('pressure (hPa)')
ax(2)=subplot(3,1,2);
pcolor(mday6h([1:end end])-date0-0.25,plev,ugrads([1:end end],[1:end end])); shading flat; axis ij
title('{\bf V}\cdot\nabla{s} (J kg^{-1} s^{-1})','fontweight','normal')
caxis([-2 2]*1.e-3);
cb = colorbar; set(cb,'tickdir','out')
ax(3)=subplot(3,1,3);
pcolor(mday6h([1:end end])-date0-0.25,plev,omgdsdp([1:end end],[1:end end])); shading flat; axis ij
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
pcolor(mday6h([1:end end])-date0-0.25,plev,Ldqbardt_ext([1:end end],[1:end end])); shading flat; axis ij
title('-L\partial{q}/\partial{t} (J kg^{-1} s^{-1})','fontweight','normal')
caxis([-2 2]*1e-1);
cb = colorbar; set(cb,'tickdir','out')
ylabel('pressure (hPa)')
ax(2)=subplot(3,1,2);
pcolor(mday6h([1:end end])-date0-0.25,plev,Lugradq([1:end end],[1:end end])); shading flat; axis ij
title('-L{\bf{V}}\cdot\nabla{q} (J kg^{-1} s^{-1})','fontweight','normal')
caxis([-2 2]*1.e-2);
cb = colorbar; set(cb,'tickdir','out')
ax(3)=subplot(3,1,3);
pcolor(mday6h([1:end end])-date0-0.25,plev,Lomgdqdp([1:end end],[1:end end])); shading flat; axis ij
title('-L\omega\partial{q}/\partial{p} (J kg^{-1} s^{-1})','fontweight','normal')
xlabel('July day')
caxis([-4 4.001]*1.e-1);
cb = colorbar; set(cb,'tickdir','out')
set(ax(:),'fontsize',15,'xlim',xlm,'ylim',ylm,'xtick',xtk,...
    'xticklabel',xtl,'ytick',ytl,'tickdir','out')
orient tall
saveas(gcf,'Q2adv_terms.png')

clf
subplot(2,2,1)
[Xl,Yl] = m_ll2xy([70 100], [0 25]);
m_proj('Equidistant Cylindrical', 'lon',[70 100],'lat',[0 25])
axis tight; axis equal;
m_gshhs_h('patch',[.8 .8 .8],'edgecolor','k'); % m_coast()
hold on
m_plot(latlon(:,2),latlon(:,1),'.','color','r','markersize',12)
%m_grid('box','fancy','tickdir','in', 'lon',[70 100],'lat',[0 25]) % fails
% work around 
xt = 70:10:100; yt = 0:10:25;
[Xt,~] = m_ll2xy(xt, 0*xt, 'patch'); [~,Yt] = m_ll2xy(0*yt, yt, 'patch');
set(gca, 'xtick', Xt, 'xticklabel', xt, 'ytick', Yt, 'yticklabel', yt)
set(gca,'fontsize',12) %, 'dataaspectratio',[1 1 1],'xlim',Xl,'ylim',Yl)
saveas(gcf,'sounding_array_plan.eps')
