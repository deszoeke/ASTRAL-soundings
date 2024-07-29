% Merge U. Wyoming formatted soundings with IMD soundings for MISO-BOB.

griddir = '/Volumes/cruise/SR1911/share/data/radiosonde/gridded/';
stations = {'Chennai','Kolkata','PortBlair','Thiruvananthapuram','Vishakhapatnam','SallyRide'};

for m = 1:length(stations)
%     fname = dir(fullfile( griddir, [ char(stations(m)) '.mat'] ));
%     wname = dir(fullfile( griddir, [ char(stations(m)) '_wyoming.mat'] ));
    snd  = getfield( load( fullfile( griddir, [ char(stations(m)) '_grd.mat'] ) ), 'grd');
    wsnd = getfield( load( fullfile( griddir, [ char(stations(m)) '_grd_wyoming.mat'] ) ), 'grd');
    
