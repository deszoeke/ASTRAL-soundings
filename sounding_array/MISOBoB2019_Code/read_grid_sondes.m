% read gridded sondes and plot time height series
addpath('/Volumes/cruise/SR1911/share/scripts/radiosonde/subs');

% read and plot gridded sounding data
anom = @(x) x-nanmean(x,2);

% load stations
datadir = '/Volumes/cruise/SR1911/share/data/radiosonde/gridded/';

stations_imd = {'Chennai','Kolkata','PortBlair','Thiruvananthapuram','Vishakhapatnam','SallyRide'};
stations_wyo = {'Bhubaneswar','Chennai','Kolkata','PortBlair','Vishakhapatnam'}; %,'Thiruvananthapuram'};
stations = union(stations_imd, stations_imd);
% stations = {'Chennai', 'Kolkata', 'PortBlair', 'SallyRide',
% 'Thiruvananthapuram', 'Vishakhapatnam' };

% Chennai: 13.0827° N, 80.2707° E
% Kolkata: 22.5726° N, 88.3639° E
% Port Blair: 11.6234° N, 92.7265° E
% Bhubaneshwar: 20.2961° N, 85.8245° E
% Sally Ride: 15.5° N, 89° E
% Thiruvananthapuram: 8.5241° N, 76.9366° E
% Vishakhapatnam: 17.6868° N, 83.2185° E
latlon = [ 13.0827, 80.2707;
           22.5726, 88.3639;
           11.6234, 92.7265;
           20.2961, 85.8245;
              15.5, 89;
            8.5241, 76.9366;
           17.6868, 83.2185; ];

% load soundings from 'Chennai','Kolkata','PortBlair', and 'Sally Ride'
for ist = 1:4
    arr(ist) = load( [fullfile( datadir, stations{ist} ), '_grd.mat'] ); % --> grd
%     grd(ist).timestr = datestr(grd(ist).time, 'yyyy-mm-dd HH:MM'); % fill in if it's missing
end
arr(1).grd.station = 'Chennai';
arr(2).grd.station = 'Kolkata';
arr(3).grd.station = 'Port Blair';
arr(4).grd.station = 'Sally Ride';

b2rcolormap(13);
[axth, cbth] = timepress_anom_sonde(arr, 'ptmp', 'potential temperature anomaly (K)', 300:5:600);
set(axth(:),'ylim', [100 1000], 'ytick',200:200:1000)
% strong diurnal cycles in Chennai, near surface but also aloft.
print('-dpng','timepress_th.png')

[axrh, cbrh] = timepress_anom_sonde(arr, 'rh', 'relative humidity anomaly (%)', 0:20:100);
set(axrh(:), 'clim', [-40 40])
print('-dpng','timepress_rh.png')

[axu, cbu] = timepress_anom_sonde(arr, 'u', 'zonal wind anomaly (m/s)', -50:10:50);
set(axu(:), 'clim', [-20 20])
print('-dpng','timepress_u.png')

[axv, cbv] = timepress_anom_sonde(arr, 'v', 'meridional wind anomaly (m/s)', -50:10:50);
set(axv(:), 'clim', [-15 15])
print('-dpng','timepress_v.png')
