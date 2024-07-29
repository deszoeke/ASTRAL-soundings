% read all Yap and Palau soundings and average them between pressure levels
clear

addpath(genpath('/Users/sdeszoek/Documents/MATLAB/user_tools/stat/'));
addpath(genpath('/Users/sdeszoek/Documents/MATLAB/user_tools/thermo/'));

% pressure levels
pint = [0:25:150 200:50:1000]';
pmid = pint(1:end-1) + diff(pint)/2;
np = length(pmid);

% load and average Yap sondes to pressure levels
Snd   = load('/Users/sdeszoek/Data/cruises/PISTON_MISOBOB_2018/Soundings/Yap/Yap_soundings.mat');
nsonde = size(Snd.hght,1);

[hght, temp, rhum, uwnd, vwnd, mixr, shum, thta, thte, dwpt] = deal(NaN(np, nsonde));
[time, lat, lon, psfc] = deal(NaN(1,nsonde));
psfc = max(Snd.pres,[],2);
time = Snd.datenumber;

cut = @(x,i) x(i);

for i = 1:nsonde
    % time(i) = nanmean(sonde.Time); % seconds since Jan 1, 2018 = (yday-1)*86400
    %[up1, uorder1] = unique(Snd.pres(i,:)');
    %up = up1(end:-1:1); uorder = uorder1(end:-1:1);
    up = Snd.pres(i,:)';
    hght(:,i) = cut(binavg(pint, up, Snd.hght(i,:)',1),1:np);
    temp(:,i) = cut(binavg(pint, up, Snd.temp(i,:)',1),1:np); % K
    rhum(:,i) = cut(binavg(pint, up, Snd.relh(i,:)',1),1:np);
    uwnd(:,i) = cut(binavg(pint, up, Snd.uwnd(i,:)',1),1:np);
    vwnd(:,i) = cut(binavg(pint, up, Snd.vwnd(i,:)',1),1:np);
    mixr(:,i) = cut(binavg(pint, up, Snd.mixr(i,:)',1),1:np);
    shum(:,i) = cut(binavg(pint, up, Snd.shum(i,:)',1),1:np);
    thta(:,i) = cut(binavg(pint, up, Snd.thta(i,:)',1),1:np);
    thte(:,i) = cut(binavg(pint, up, Snd.thte(i,:)',1),1:np);
    dwpt(:,i) = cut(binavg(pint, up, Snd.dwpt(i,:)',1),1:np);
end

% Pi = temp./thta;
thdw = dwpt.*(thta./temp); % dew point potential temperature
thes = theta_es(pmid*100, temp); % saturation equivalent potential temperature
qast = 1.0e3*qs(pmid*100,temp-273.15); % saturation specific humidity
thtw = theta_w(pmid*100, temp, shum/1e3); % interesting how constant over troposphere
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

temp = temp + 273.15;

clear Snd cut up order i
save('Yap_sondes_pavg.mat',...
     'time','yday','lat' ,'lon' ,'psfc','pint',...
     'hght','temp','rhum','uwnd','vwnd',...
     'shum','thta','thte','dwpt','thdw','thes','qast','thtw'    );

% load and average Koror sondes to pressure levels
Snd = load('/Users/sdeszoek/Data/cruises/PISTON_MISOBOB_2018/Soundings/Palau/91408_KororPalau/Palau_soundings.mat');
nsonde = size(Snd.hght,1);

[hght, temp, rhum, uwnd, vwnd, shum, thta, thte, dwpt] = deal(NaN(np, nsonde));
[time, lat, lon, psfc] = deal(NaN(1,nsonde));
psfc = max(Snd.pres,[],2);
time = Snd.datenumber;

cut = @(x,i) x(i);

for i = 1:nsonde    
%     [up, uorder] = unique(Snd.pres(:,i));
    up = Snd.pres(i,:)';
    hght(:,i) = cut(binavg(pint, up, Snd.hght(i,:)'),1:np);
    temp(:,i) = cut(binavg(pint, up, Snd.temp(i,:)'),1:np); % K
    rhum(:,i) = cut(binavg(pint, up, Snd.relh(i,:)'),1:np);
    uwnd(:,i) = cut(binavg(pint, up, Snd.uwnd(i,:)'),1:np);
    vwnd(:,i) = cut(binavg(pint, up, Snd.vwnd(i,:)'),1:np);
    mixr(:,i) = cut(binavg(pint, up, Snd.mixr(i,:)'),1:np);
    shum(:,i) = cut(binavg(pint, up, Snd.shum(i,:)'),1:np);
    thta(:,i) = cut(binavg(pint, up, Snd.thta(i,:)'),1:np);
    thte(:,i) = cut(binavg(pint, up, Snd.thte(i,:)'),1:np);
    dwpt(:,i) = cut(binavg(pint, up, Snd.dwpt(i,:)'),1:np);
end

% Pi = temp./thta;
thdw = dwpt.*(thta./temp); % dew point potential temperature
thes = theta_es(pmid*100, temp); % saturation equivalent potential temperature
qast = 1.0e3*qs(pmid*100,temp-273.15); % saturation specific humidity
thtw = theta_w(pmid*100, temp, shum/1e3); % interesting how constant over troposphere
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

temp = temp + 273.15;

% clear Snd cut up uorder i
save('Palau_sondes_pavg.mat',...
     'time','yday','lat' ,'lon' ,'yday','psfc','pint',...
     'hght','temp','rhum','uwnd','vwnd',...
     'shum','thta','thte','dwpt','thes','qast','thtw'    );

