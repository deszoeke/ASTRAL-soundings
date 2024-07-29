% compare integrated vapor WVP from microwave with soundings

filepath = '/Users/sdeszoek/Data/cruises/PISTON_MISOBOB_2018/Soundings/Thompson/netcdf';
flist = dir(fullfile(filepath,'Thompson*.nc'));
nsonde = length(flist);

%timestamp = sprintf('%4d%02d%02d_%02d%02d%02d', year, month, day, hour, minute, second);
%ncf = fullfile(filepath, ['Thompson_sonde' timestamp '.nc']);

zint = (0:10:27000)';
nz = length(zint);
[pres, temp, rhum, uwnd, vwnd, shum, thta, thte, dwpt, wspd, wdir] = deal(NaN(nz, nsonde));
[time, lat, lon] = deal(NaN(1,nsonde));
for i = 1:nsonde
    ncf = fullfile(filepath, flist(i).name);
    
    sonde = read_TGT_sonde(ncf);
    time(i) = nanmean(sonde.Time); % seconds since Jan 1, 2018 = (yday-1)*86400
    lat(i) = nanmean(sonde.Latitude_of_balloon); % seconds since Jan 1, 2018 = (yday-1)*86400
    lon(i) = nanmean(sonde.Longitude_of_balloon); % seconds since Jan 1, 2018 = (yday-1)*86400
    [t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11] = interp_sonde_z(sonde, zint);
    pres(:,i) = t1;
    temp(:,i) = t2; % K
    rhum(:,i) = t3;
    uwnd(:,i) = t4;
    vwnd(:,i) = t5;
    shum(:,i) = t6;
    thta(:,i) = t7;
    thte(:,i) = t8;
    dwpt(:,i) = t9;
    wspd(:,i) = t10;
    wdir(:,i) = t11;

end

% Pi = temp./thta;
thdw = dwpt.*(thta./temp); % dew point potential temperature
thes = theta_es(pres*100, temp); % saturation equivalent potential temperature
qast = 1.0e3*qs(pres*100,temp-273.15); % saturation specific humidity

thtw = theta_w(pres*100, temp, shum/1e3); % interesting how constant over troposphere

yday = time/86400 + 1;

save('TGT_sondes_z.mat',...
     'time','lat' ,'lon' ,'yday','zint',...
     'pres','temp','rhum','uwnd','vwnd',...
     'shum','thta','thte','dwpt','wspd',...
     'wdir','thdw','thes','qast','thtw'    );
