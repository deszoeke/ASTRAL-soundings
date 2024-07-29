% read all TGT soundings and average them between pressure levels
clear all;

addpath(genpath('/Users/sdeszoek/Documents/MATLAB/user_tools/stat/'));
addpath(genpath('/Users/sdeszoek/Documents/MATLAB/user_tools/thermo/'));

filepath = '/Users/sdeszoek/Data/cruises/PISTON_MISOBOB_2018/Soundings/Thompson/netcdf';
flist = dir(fullfile(filepath,'Thompson*.nc'));
nsonde = length(flist);

%timestamp = sprintf('%4d%02d%02d_%02d%02d%02d', year, month, day, hour, minute, second);
%ncf = fullfile(filepath, ['Thompson_sonde' timestamp '.nc']);
%thermo_constants;

pint = [0:25:150 200:50:1000]';
pmid = pint(1:end-1) + diff(pint)/2;
np = length(pmid);

[hght, temp, rhum, uwnd, vwnd, shum, thta, thte, dwpt] = deal(NaN(np, nsonde));
[time, lat, lon, psfc] = deal(NaN(1,nsonde));
for i = 1:nsonde
    ncf = fullfile(filepath, flist(i).name);
    
    sonde = read_TGT_sonde(ncf);
    time(i) = nanmean(sonde.Time); % seconds since Jan 1, 2018 = (yday-1)*86400
    lat(i)  = nanmean(sonde.Latitude_of_balloon); % seconds since Jan 1, 2018 = (yday-1)*86400
    lon(i)  = nanmean(sonde.Longitude_of_balloon); % seconds since Jan 1, 2018 = (yday-1)*86400
    psfc(i) = sonde.Release_Pressure;
    [t1, t2, t3, t4, t5, t6, t7, t8, t9] = avg_sonde_p(sonde, pint);
    hght(:,i) = t1(1:end-1);
    temp(:,i) = t2(1:end-1); % K
    dwpt(:,i) = t3(1:end-1);
    rhum(:,i) = t4(1:end-1);
    shum(:,i) = t5(1:end-1)*1.e-3;
    thta(:,i) = t6(1:end-1);
    thte(:,i) = t7(1:end-1);
    uwnd(:,i) = t8(1:end-1);
    vwnd(:,i) = t9(1:end-1);
end

% Pi = temp./thta;
thdw = dwpt.*(thta./temp); % dew point potential temperature
thes = theta_es(pmid*100, temp); % saturation equivalent potential temperature
qast = qs(pmid*100,temp-273.15); % saturation specific humidity

thtw = theta_w(repmat(pmid,[1 nsonde])*100, temp, shum); % interesting how constant over troposphere

yday = time/86400 + 1;
% For some reason 20180906_173022 (249.7687) is a duplicate of 20180906_143022.
% A look at the edt files shows no 17:30:22 sounding.
% 20180906_173022 is removed by renaming it to 
% duplicates_20180906_143022_Thompson_sonde_20180906_173022.nc
% and the launch time for yday(126) should be 14:30:22 (249.6044).
if yday(126) > 249.75 % hack until fixed in netcdf files
    yday(126) = datenum(2018,09,06,14,30,22)-datenum(2018,0,0);
end
% time is still wrong?

clear t1 t2 t3 t4 t5 t6 t7 t8 t9 filepath flist i ncf
save('TGT_sondes_pavg.mat',...
     'time','yday','lat' ,'lon' ,'psfc','pint',...
     'hght','temp','rhum','uwnd','vwnd',...
     'shum','thta','thte','dwpt','thdw','thes','qast','thtw'    );
