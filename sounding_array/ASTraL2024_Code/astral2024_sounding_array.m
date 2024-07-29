% Calculate divergence, vorticity, and omega from sounding array for ASTRAL IOP1
% Alex Kinsella, June 2024

%% Define stations and their lat/lons
Stn = ["ahmedabad","bhubaneswar","chennai","karaikal","kolkata","kochi","mangalore","minicoy","portblair","pune","santacruz","tgt","visakhapatnam"];
nstn = numel(Stn);
itgt = find(Stn=="tgt"); % Station number for TGT

pres = 1020:-5:0; % Pressure levels that files are binned to 
np = numel(pres);

tvec = datetime(2024,4,29,0,0,0):hours(2):datetime(2024,6,4,0,0,0); % 2 hour time grid for reading in soundings
nt = numel(tvec);

load("/Volumes/cruiseshare/For_Science/Situational_Awareness_ShipboardData/gps_compiled.mat");
tgtdt = datetime(2024,1,1) + days(GPS_c.dday);
tgtlat = interp1(tgtdt,GPS_c.lat,tvec); % Get TGT locations on time grid 
tgtlon = interp1(tgtdt,GPS_c.lon,tvec);

latlon = nan(nstn,2,nt);
latlon(1,1,:) = 23.0225; latlon(1,2,:) = 72.5714; % Ahmedabad
latlon(2,1,:) = 20.2961; latlon(2,2,:) = 85.8245; % Bhubaneswar
latlon(3,1,:) = 13.0827; latlon(3,2,:) = 80.2707; % Chennai
latlon(4,1,:) = 10.9254; latlon(4,2,:) = 79.8380; % Karaikal
latlon(5,1,:) = 9.9312; latlon(5,2,:) = 76.2673; % Kochi
latlon(6,1,:) = 22.5726; latlon(6,2,:) = 88.3639; % Kolkata
latlon(7,1,:) = 12.9141; latlon(7,2,:) = 74.8560; % Mangalore 
latlon(8,1,:) = 8.2740; latlon(8,2,:) = 73.0496; % Minicoy
latlon(9,1,:) = 11.6234; latlon(9,2,:) = 92.7265; % Port Blair
latlon(10,1,:) = 18.5204; latlon(10,2,:) = 73.8567; % Pune
latlon(11,1,:) = 19.0843; latlon(11,2,:) = 72.8360; % Santa Cruz
latlon(12,1,:) = tgtlat; latlon(12,2,:) = tgtlon; % TGT
latlon(13,1,:) = 17.6868; latlon(13,2,:) = 83.2185; % Vishakhapatnam

[u_all,v_all,T_all,z_all,q_all] = deal(nan(np,nt,nstn));
% Read in files
%files = dir(["/Users/kinsella/Documents/WHOI/Projects/ASTraL/2024_IOP1/ASTRAL-soundings/data/grid_sondes/*.nc"]);
files = dir(["/Volumes/cruiseshare/For_Science/Situational_Awareness_Processing/data/radiosonde/grid_sondes/*.nc"]);
for ff = 1:numel(files)
    if strcmp(files(ff).name(1:3),'TGT')
        time = datetime(str2num(files(ff).name(5:8)),str2num(files(ff).name(9:10)),str2num(files(ff).name(11:12)),str2num(files(ff).name(13:14)),0,0);
        itime = find(time==tvec);
        if ~isempty(itime)
            u_all(:,itime,itgt) = ncread([files(ff).folder,'/',files(ff).name],'u'); 
            v_all(:,itime,itgt) = ncread([files(ff).folder,'/',files(ff).name],'v'); 
            T_all(:,itime,itgt) = ncread([files(ff).folder,'/',files(ff).name],'temp'); 
            z_all(:,itime,itgt) = ncread([files(ff).folder,'/',files(ff).name],'hght'); 
            if time<datetime(2024,5,16) % q was in kg/kg for leg 1 and g/kg for leg 2. Update if this gets standardized. 
                q_all(:,itime,itgt) = 1000*ncread([files(ff).folder,'/',files(ff).name],'q'); 
            else
                q_all(:,itime,itgt) = ncread([files(ff).folder,'/',files(ff).name],'q'); 
            end
        end
    else
        time = datetime(str2num(files(ff).name(1:4)),str2num(files(ff).name(5:6)),str2num(files(ff).name(7:8)),str2num(files(ff).name(10:11)),0,0);
        name = lower(files(ff).name(13:end-3));
        istn = find(Stn==name);
        itime = find(time==tvec);
        if ~isempty(itime)
            u_all(:,itime,istn) = ncread([files(ff).folder,'/',files(ff).name],'u'); 
            v_all(:,itime,istn) = ncread([files(ff).folder,'/',files(ff).name],'v'); 
            T_all(:,itime,istn) = ncread([files(ff).folder,'/',files(ff).name],'T'); 
            z_all(:,itime,istn) = ncread([files(ff).folder,'/',files(ff).name],'z'); 
            Td = ncread([files(ff).folder,'/',files(ff).name],'Td'); % Dewpoint
            es = 6.112*exp((17.67*Td)./(Td+243.5)); % Saturation vapor pressure 
            w = 0.622*es./((pres'-es)); % Mixing ratio 
            q_all(:,itime,istn) = 1000*w./(1+w); % Specific humidity in g/kg
        end
    end
end

%% Fix bad TGT sounding
ibad = find(tvec == datetime(2024,5,8,0,0,0));
prescheck = pres<850;
u_all(prescheck,ibad,itgt) = nan;
v_all(prescheck,ibad,itgt) = nan;
T_all(prescheck,ibad,itgt) = nan;
z_all(prescheck,ibad,itgt) = nan;
q_all(prescheck,ibad,itgt) = nan;

%% Bin to coarser levels
dp = 50; % Bin width (mb)
pbins = 1000:-dp:0;
npbin = numel(pbins);
[u_bin,v_bin,T_bin,z_bin,q_bin] = deal(nan(npbin,nt,nstn));
for ip = 1:numel(pbins)
    pcheck = pres>=pbins(ip)-dp/2 & pres<pbins(ip)+dp/2;
    u_bin(ip,:,:) = nanmean(u_all(pcheck,:,:),1);
    v_bin(ip,:,:) = nanmean(v_all(pcheck,:,:),1);
    T_bin(ip,:,:) = nanmean(T_all(pcheck,:,:),1);
    z_bin(ip,:,:) = nanmean(z_all(pcheck,:,:),1);
    q_bin(ip,:,:) = nanmean(q_all(pcheck,:,:),1);
end

%% Interpolate to 6-hourly 
t_grid = datetime(2024,4,29,0,0,0):hours(6):datetime(2024,6,4,0,0,0);
ntgrid = numel(t_grid);

[u_interp,v_interp,T_interp,z_interp,q_interp] = deal(nan(npbin,ntgrid,nstn));
for istn = 1:nstn
    for ip = 1:npbin
        goodcheck = ~isnan(nanmean(u_bin(:,:,istn),1));
        u_interp(ip,:,istn) = interp1(tvec(goodcheck),squeeze(u_bin(ip,goodcheck,istn)),t_grid);
        %u_fill(ip,:,istn) = fillmissing(u_interp(ip,:,istn),'linear');
        v_interp(ip,:,istn) = interp1(tvec(goodcheck),squeeze(v_bin(ip,goodcheck,istn)),t_grid);
        %v_fill(ip,:,istn) = fillmissing(v_interp(ip,:,istn),'linear');
        T_interp(ip,:,istn) = interp1(tvec(goodcheck),squeeze(T_bin(ip,goodcheck,istn)),t_grid);
        z_interp(ip,:,istn) = interp1(tvec(goodcheck),squeeze(z_bin(ip,goodcheck,istn)),t_grid);
        q_interp(ip,:,istn) = interp1(tvec(goodcheck),squeeze(q_bin(ip,goodcheck,istn)),t_grid);
    end
end

% Remove interpolated values that have no soundings within a certain number of hours
% Do independently for each variable in case of sensor failures
nhrs = 24; % Tolerance for extended interpolation (hrs)
for istn = 1:nstn
    for it = 1:ntgrid
        tcheck = abs(tvec-t_grid(it)) < hours(nhrs);
        good = ~isnan(nanmean(u_bin(:,tcheck,istn),'all'));
        if good == false
            u_interp(:,it,istn) = nan;
        end

        good = ~isnan(nanmean(v_bin(:,tcheck,istn),'all'));
        if good == false
            v_interp(:,it,istn) = nan;
        end

        good = ~isnan(nanmean(T_bin(:,tcheck,istn),'all'));
        if good == false
            T_interp(:,it,istn) = nan;
        end

        good = ~isnan(nanmean(z_bin(:,tcheck,istn),'all'));
        if good == false
            z_interp(:,it,istn) = nan;
        end

        good = ~isnan(nanmean(q_bin(:,tcheck,istn),'all'));
        if good == false
            q_interp(:,it,istn) = nan;
        end
    end
end 


%% Make windspeed plot for each station
fig1 = makefig;
for istn = 1:nstn
    subplot(4,4,istn)
    pcolor(t_grid,pbins,sqrt(u_interp(:,:,istn).^2+v_interp(:,:,istn).^2))
    shading flat
    set(gca,'ydir','reverse')
    makespruce(20)
    maketitle(Stn(istn),20)
end

%% Make specific humidity plot for each station
fig1 = makefig;
for istn = 1:nstn
    subplot(4,4,istn)
    pcolor(t_grid,pbins,q_interp(:,:,istn))
    cb = makecb([0 30],cmocean('rain'),'Specific Humidity (g/kg)',30);
    shading flat
    set(gca,'ydir','reverse')
end

%% Make temperature plot for each station
fig1 = makefig;
for istn = 1:nstn
    subplot(4,4,istn)
    pcolor(t_grid,pbins,T_interp(:,:,istn))
    cb = makecb([-60 30],cmocean('thermal'),'Temperature (deg C)',30);
    shading flat
    set(gca,'ydir','reverse')
end


%% Calculate divergence 
% choose a reference latitude
lat = squeeze(latlon(:,1,:));
lon = squeeze(latlon(:,2,:));
coslat = cosd(lat); 
sinlat = sind(lat);
r_earth = 6.371e6; % m

y = r_earth * lat;          % meridional coordinate, meters
x = r_earth * coslat .* lon; % zonal coordinate, meters

% Fit plane to measurements to calculate derivatives 
[dudx,dudy,dvdx,dvdy] = deal(nan(npbin,ntgrid));
for it = 1:ntgrid
    for ip = 1:npbin
        dudx(ip,it)  = partiald(x(:,it)', squeeze( u_interp(ip,it,:))', 2)';
        dudy(ip,it)  = partiald(y(:,it)', squeeze( u_interp(ip,it,:))', 2)';
        dvdx(ip,it)  = partiald(x(:,it)', squeeze( v_interp(ip,it,:))', 2)';
        dvdy(ip,it)  = partiald(y(:,it)', squeeze( v_interp(ip,it,:))', 2)';
    end
end

divg = dudx + dvdy; % Divergence 
vort = dvdx - dudy; % Vorticity

%% Plot divergence and vorticity
fig1 = makefig;
pcolor(t_grid,pbins,divg)
cb = makecb([-2e-7 2e-7],cmocean('balance'),'Divergence (1/s)',30);
set(gca,'ydir','reverse')
shading flat
makespruce(30)
maketitle('Divergence During EKAMSAT IOP1',40)

fig2 = makefig;
pcolor(t_grid,pbins,vort)
cb = makecb([-2e-7 2e-7],cmocean('balance'),'Vorticity (1/s)',30);
set(gca,'ydir','reverse')
shading flat
makespruce(30)
maketitle('Vorticity During EKAMSAT IOP1',40)

%% Calculate and plot omega 
% adjust mean column divergence to be 0 to ensure 
% column-integrated nondivergence and omega(sfc) = omega(p=0) = 0
divg_integrand = divg - nanmean(divg, 1);
% Integrate omega up from zero at the surface
% Pa/s on pressure level interfaces, starting at pint(1)= 0 hPa
%omga = flipud([zeros(1,ntgrid); cumsum(flipud( 100*dp .* divg_integrand ))]);
omga = flipud([cumsum(flipud( 100*(-dp) .* divg_integrand ),'omitnan')]);

fig3 = makefig;
pcolor(t_grid,pbins,omga)
cb = makecb([-5e-3 5e-3],cmocean('balance'),'Pressure Velocity (Pa/s)',30);
set(gca,'ydir','reverse')
shading flat
makespruce(30)
maketitle('Large-Scale Omega During EKAMSAT IOP1',40)


%% Plot div,vort,omega

fig1 = makefig;

ax1 = subplot(2,1,1);
pcolor(t_grid,pbins,divg)
cb = makecb([-3e-7 3e-7],cmocean('balance'),'Divergence (1/s)',20,ax1);
set(gca,'ydir','reverse')
shading flat
makespruce(30)
maketitle('Divergence During EKAMSAT IOP1',40)

ax2 = subplot(2,1,2);
pcolor(t_grid,pbins,vort)
cb = makecb([-3e-7 3e-7],cmocean('balance'),'Vorticity (1/s)',20,ax2);
set(gca,'ydir','reverse')
shading flat
makespruce(30)
maketitle('Vorticity During EKAMSAT IOP1',40)

fig2 = makefig;
pcolor(t_grid,pbins,omga)
cb = makecb([-5e-3 5e-3],cmocean('balance'),'Pressure Velocity (Pa/s)',30);
set(gca,'ydir','reverse')
shading flat
makespruce(30)
maketitle('Large-Scale Omega During EKAMSAT IOP1',40)


%% Make plot of all stations

load coastlines
fig1 = makefig;
scatter(latlon(:,2,1),latlon(:,1,1),50,'filled','color','r')
hold on
plot(coastlon,coastlat)
xlim([70 100])
ylim([0 30])
daspect([1 cos(pi/180*10) 1])