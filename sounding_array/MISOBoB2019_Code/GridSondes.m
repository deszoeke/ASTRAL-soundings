%% Grid Radiosondes %%
% First run ReadSoundings, ReadWyoSoundings then:
% GridSondes merges U. Wyoming, IMD, and Sally Ride soundings into a single
% file per station.
clear all; close all;

%% set paths
sysType = computer;
username = char(java.lang.System.getProperty('user.name'));
if strncmp(sysType,'MACI64',7)     % set Mac paths
    switch username
        case 'sdeszoek'
            datadir = '/Users/sdeszoek/Data/cruises/MISOBOB_2019/SR1911/radiosonde/';
            addpath('~/Data/cruises/MISOBOB_2019/SR1911/scientific_analysis/programs/radiosonde/subs/');
            addpath('~/Data/cruises/MISOBOB_2019/SR1911/scientific_analysis/programs/radiosonde/radiosonde_PISTON2018/thermo/')
        otherwise
            datadir = '/Volumes/cruise/SR1911/share/data/radiosonde/';
            addpath('/Volumes/cruise/SR1911/share/scripts/radiosonde/subs/');
            addpath('/Volumes/cruise/SR1911/share/scripts/radiosonde/radiosonde_PISTON2018/thermo/')
    end
end
matdir = fullfile(datadir, '/mat/');

% imd_stations={'Chennai','Kolkata','PortBlair','Thiruvananthapuram','Vishakhapatnam'};
stations_imd = {'Chennai','Kolkata','PortBlair','Thiruvananthapuram','Vishakhapatnam','SallyRide'};
stations_wyo = {'Bhubaneswar','Chennai','Kolkata','PortBlair','Vishakhapatnam'}; %,'Thiruvananthapuram'};
% stations in _both_ Wyo. and IMD
stations_mrg = intersect(stations_imd, stations_wyo);

% stations that are _only_ in Wyo. xor IMD collections
[stns_xor, io_imd, io_wyo] = setxor(stations_imd, stations_wyo);
stns_xor_imd = stations_imd(io_imd);
stns_xor_wyo = stations_wyo(io_wyo);

dp = 5; % set pressure resolution
p = (1015:-dp:10)';
np = length(p);
trunc = @(x) flipud(x(1:np));

% constants
Rd = 287.04;
Cp = 1005.7;
kts2ms = 0.51444444; % all snd files are already converted to m/s

%% grid stations with IMD and U. Wyo
for istn = 1:length(stations_mrg)
    
    fnames = dir([matdir stations_mrg{istn} '_*.mat']);
    fnames_wyo = dir( [fullfile(matdir, 'wyoming', stations_mrg{istn}) '_*.mat'] );
    
    tmp = cell2mat({fnames.name}');
    time_string = tmp(:,end-15:end-4);
    time = datenum(time_string, 'yyyymmddHHMM');
    tmp = cell2mat({fnames_wyo.name}');
    time_string_wyo = tmp(:,end-15:end-4);
    time_wyo = datenum(time_string_wyo, 'yyyymmddHHMM');
    
    % add (many) Wyo soundings not already in IMD soundings
    % logical selects Wyo soundings
    % find unique times separated by >2 hours
    iwyo = all ( abs(time-time_wyo') > 2/24 );
    nwyo = sum(iwyo);
    fwyo = find( iwyo );
    
    % concatenate time
    time_mrg1 = [time(:); time_wyo(iwyo)];
    [time_mrg, order] = sort(time_mrg1);
    
    % initialize grd structure
    grd   = [];
    grd.p = p;
    grd.time = nan(1,  length(fnames));
%     grd.timestr = [];
    grd.source = {};
    grd.temp = nan(np, length(fnames));
    grd.ptmp = nan(np, length(fnames));
    grd.etmp = nan(np, length(fnames));
    grd.thw  = nan(np, length(fnames));
    grd.dew  = nan(np, length(fnames));
    grd.rh   = nan(np, length(fnames));
    grd.qv   = nan(np, length(fnames));
    grd.h    = nan(np, length(fnames));
    grd.gp   = nan(np, length(fnames));
    grd.u    = nan(np, length(fnames));
    grd.v    = nan(np, length(fnames));
    grd.lat  = nan(np, length(fnames));
    grd.lon  = nan(np, length(fnames));
    
    % load IMD soundings
    for ii = 1:length(fnames)
        load( fullfile( fnames(ii).folder, fnames(ii).name ) ); % snd
%         wyo = getfield(load( fullfile( fnames(ii).folder, fnames(ii).name ) ), 'snd');
        %% derived quantities
        qsat = 1.0e3 * qs(snd.p*100.0, snd.t);
        
        % sonde.Specific_Humidity = 0.01*sonde.Relative_Humidity.*qsat([sonde.Temperature-273.15, sonde.Pressure]);
        if ~isfield(snd, 'qv')
            snd.qv = 0.01 * snd.rh .* qsat;
        end
        
        Pi = (snd.p*1e-3).^(Rd/Cp*(1-0.28*snd.qv*1.0e-3));
        snd.ptmp = (snd.t+273.15)./Pi;
        
        Tl = Tlcl(snd.p*1.0e2, snd.t+273.15, snd.qv);
        
        snd.etmp = snd.ptmp .* exp( (3376./Tl - 2.54) .* ...
            snd.qv*1.0e-3 .* (1 + 0.81*snd.qv*1.0e-3) ); % Bolton 1980
        
        snd.u = -snd.ws.*sind(double(snd.wd));
        snd.v = -snd.ws.*cosd(double(snd.wd));
        
        snd.thw = theta_w(snd.p*100,snd.t+273.15,snd.qv./1000);
        
        grd.time(ii) = time(ii); %snd.time;
        disp( datestr( time(find(ii, 1, 'last')) ) );
        
        grd.source(ii) = {snd.source};
        
        % average soundings to dp intervals where possible
        pbin = [10-dp/2; dp/2+flipud(grd.p)];
        grd.temp(:,ii) = trunc(binavg( pbin, snd.p, snd.t)   );
        grd.ptmp(:,ii) = trunc(binavg( pbin, snd.p, snd.ptmp));
        grd.etmp(:,ii) = trunc(binavg( pbin, snd.p, snd.etmp));
        grd.thw( :,ii) = trunc(binavg( pbin, snd.p, snd.thw) );
        grd.dew( :,ii) = trunc(binavg( pbin, snd.p, snd.td)  );
        grd.qv(  :,ii) = trunc(binavg( pbin, snd.p, snd.qv)  );
        grd.u(   :,ii) = trunc(binavg( pbin, snd.p, snd.u)   );
        grd.v(   :,ii) = trunc(binavg( pbin, snd.p, snd.v)   );
        grd.rh(  :,ii) = trunc(binavg( pbin, snd.p, snd.rh)  );
        if isfield(snd,'h')
            grd.h(  :,ii) = trunc(binavg( pbin, snd.p, snd.h)   );
        else
            grd.h(  :,ii) = NaN;
        end
        if isfield(snd,'gp')
            grd.gp( :,ii) = trunc(binavg( pbin, snd.p, snd.gp)  );
        else
            grd.gp(  :,ii) = NaN;
        end
        if isfield(snd,'lat')
            grd.lat(:,ii) = trunc(binavg(pbin, snd.p, snd.lat) );
            grd.lon(:,ii) = trunc(binavg(pbin, snd.p, snd.lon) );
        else
            grd.lat(:,ii) = NaN;
            grd.lon(:,ii) = NaN;
        end
        
        % otherwise interpolate to p grid
        [pun, un] = unique(snd.p);
        if abs(median(diff(pun))) > 2.0         
            isn = isnan(grd.temp(:,ii));
            grd.temp(isn,ii) = interp1(pun, snd.t(un)   , grd.p(isn), 'linear');
            grd.ptmp(isn,ii) = interp1(pun, snd.ptmp(un), grd.p(isn), 'linear');
            isn = isnan(grd.qv(:,ii));
            grd.etmp(isn,ii) = interp1(pun, snd.etmp(un), grd.p(isn), 'linear');
            grd.thw( isn,ii) = interp1(pun, snd.thw(un) , grd.p(isn), 'linear');
            grd.dew( isn,ii) = interp1(pun, snd.td(un)  , grd.p(isn), 'linear');
            grd.qv(  isn,ii) = interp1(pun, snd.qv(un)  , grd.p(isn), 'linear');
            isn = isnan(grd.u(:,ii));
            grd.u(   isn,ii) = interp1(pun, snd.u(un)   , grd.p(isn), 'linear');
            grd.v(   isn,ii) = interp1(pun, snd.v(un)   , grd.p(isn), 'linear');
            grd.rh(  isn,ii) = interp1(pun, snd.rh(un)  , grd.p(isn), 'linear');
            if isfield(snd,'h')
                grd.h(  :,ii) = trunc(binavg( pbin, snd.p, snd.h)   );
            else
                grd.h(  :,ii) = NaN;
            end
            if isfield(snd,'gp')
                grd.gp( :,ii) = trunc(binavg( pbin, snd.p, snd.gp)  );
            else
                grd.gp(  :,ii) = NaN;
            end
            if isfield(snd,'lat')
                isn = isnan(grd.lat(:,ii));
                grd.lat(isn,ii) = interp1(pun, snd.lat(un), grd.p(isn), 'linear');
                grd.lon(isn,ii) = interp1(pun, snd.lon(un), grd.p(isn), 'linear');
            end
        end
    end % ii - IMD soundings
    % mm stays fixed at number of IMD soundings
    mm=ii;
    
    % load and grid U. Wyoming soundings
    for nn = 1:nwyo
        % fwyo(nn) indexes fnames_wyo
        % add mm to nn to fill grd; grd needs to be put into time order
        ii = mm + nn;
        snd = getfield(load( fullfile( fnames_wyo(fwyo(nn)).folder, fnames_wyo(fwyo(nn)).name ) ), 'snd');
        %% derived quantities
        qsat = 1.0e3 * qs(snd.p*100.0, snd.t);
        
        % sonde.Specific_Humidity = 0.01*sonde.Relative_Humidity.*qsat([sonde.Temperature-273.15, sonde.Pressure]);
        if ~isfield(snd, 'qv')
            snd.qv = 0.01 * snd.rh .* qsat;
        end
        
        Pi = (snd.p*1e-3).^(Rd/Cp*(1-0.28*snd.qv*1.0e-3));
        snd.ptmp = (snd.t+273.15)./Pi;
        
        Tl = Tlcl(snd.p*1.0e2, snd.t+273.15, snd.qv);
        
        snd.etmp = snd.ptmp .* exp( (3376./Tl - 2.54) .* ...
            snd.qv*1.0e-3 .* (1 + 0.81*snd.qv*1.0e-3) ); % Bolton 1980
        
        snd.u = -snd.ws.*sind(double(snd.wd));
        snd.v = -snd.ws.*cosd(double(snd.wd));
        
        snd.thw = theta_w(snd.p*100,snd.t+273.15,snd.qv./1000);
        
        grd.time(ii) = snd.time;
        disp( datestr( snd.time(end) ) );
        
        grd.source(ii) = {'UWyo'}; % {snd.source}; % 'UWyo';
        
        % average soundings to dp intervals where possible
%         pbin = [10-dp/2; dp/2+flipud(grd.p)];
        grd.temp(:,ii) = trunc(binavg( pbin, snd.p, snd.t)    );
        grd.ptmp(:,ii) = trunc(binavg( pbin, snd.p, snd.ptmp) );
        grd.etmp(:,ii) = trunc(binavg( pbin, snd.p, snd.etmp) );
        grd.thw( :,ii) = trunc(binavg( pbin, snd.p, snd.thw)  );
        grd.dew( :,ii) = trunc(binavg( pbin, snd.p, snd.td)   );
        grd.qv(  :,ii) = trunc(binavg( pbin, snd.p, snd.qv)   );
        grd.u(   :,ii) = trunc(binavg( pbin, snd.p, snd.u)    );
        grd.v(   :,ii) = trunc(binavg( pbin, snd.p, snd.v)    );
        grd.rh(  :,ii) = trunc(binavg( pbin, snd.p, snd.rh)   );
        if isfield(snd,'h')
            grd.h(  :,ii) = trunc(binavg( pbin, snd.p, snd.h)   );
        else
            grd.h( :,ii) = NaN;
        end
        if isfield(snd,'gp')
            grd.gp( :,ii) = trunc(binavg( pbin, snd.p, snd.gp)  );
        else
            grd.gp( :,ii) = NaN;
        end
        if isfield(snd,'lat')
            grd.lat(:,ii) = trunc(binavg(pbin, snd.p, snd.lat)  );
            grd.lon(:,ii) = trunc(binavg(pbin, snd.p, snd.lon)  );
        else
            grd.lat(:,ii) = NaN;
            grd.lon(:,ii) = NaN;
        end
        
        % otherwise interpolate to p grid
        [pun, un] = unique(snd.p);
        if abs(median(diff(pun))) > 2.0         
            isn = isnan(grd.temp(:,ii));
            grd.temp(isn,ii) = interp1(pun, snd.t(un)   , grd.p(isn), 'linear');
            grd.ptmp(isn,ii) = interp1(pun, snd.ptmp(un), grd.p(isn), 'linear');
            isn = isnan(grd.qv(:,ii));
            grd.etmp(isn,ii) = interp1(pun, snd.etmp(un), grd.p(isn), 'linear');
            grd.thw( isn,ii) = interp1(pun, snd.thw(un) , grd.p(isn), 'linear');
            grd.dew( isn,ii) = interp1(pun, snd.td(un)  , grd.p(isn), 'linear');
            grd.qv(  isn,ii) = interp1(pun, snd.qv(un)  , grd.p(isn), 'linear');
            isn = isnan(grd.u(:,ii));
            grd.u(   isn,ii) = interp1(pun, snd.u(un)   , grd.p(isn), 'linear');
            grd.v(   isn,ii) = interp1(pun, snd.v(un)   , grd.p(isn), 'linear');
            grd.rh(  isn,ii) = interp1(pun, snd.rh(un)  , grd.p(isn), 'linear');
            if isfield(snd,'h')
                grd.h( isn,ii) = interp1( pun, snd.h(un), grd.p(isn), 'linear');
            else
                grd.h( :,ii) = NaN;
            end
            if isfield(snd,'gp')
                grd.gp(isn,ii) = interp1( pun, snd.gp(un), grd.p(isn), 'linear');
            else
                grd.gp( :,ii) = NaN;
            end
            if isfield(snd,'lat')
                isn = isnan(grd.lat(:,ii));
                grd.lat(isn,ii) = interp1(pun, snd.lat(un), grd.p(isn), 'linear');
                grd.lon(isn,ii) = interp1(pun, snd.lon(un), grd.p(isn), 'linear');
            end
        end
        
    end
    
    % reorder all the variables to be in order of time
    [x, ord] = sort(grd.time);
    var = fieldnames(grd);
    for i = 1:length(var)
        if size(grd.(var{i}),2) > 1
            grd.(var{i}) = grd.(var{i})(:,ord);
        end
    end
    
    grd.timestr = datestr(grd.time, 'yyyy-mm-dd HH:MM');
    
    save(fullfile(datadir, 'gridded', [stations_mrg{istn} '_grd']), 'grd');

end % IMD-Wyo merge stations

%% grid stations in Wyo only
for istn = 1:length(stns_xor_wyo)
    
    fnames = dir([fullfile(matdir, 'wyoming', stns_xor_wyo{istn}) '_*.mat']);
    
    tmp = cell2mat({fnames.name}');
    time_string = tmp(:,end-15:end-4);
    time = datenum(time_string, 'yyyymmddHHMM');
    
    % initialize grd structure    % spd moved out of load files loop
    grd   = [];
    grd.p = (1015:-dp:10)';
    np = length(grd.p);
    grd.time = nan(1,  length(fnames));
    grd.timestr = [];
    grd.temp = nan(np, length(fnames));
    grd.ptmp = nan(np, length(fnames));
    grd.etmp = nan(np, length(fnames));
    grd.thw  = nan(np, length(fnames));
    grd.dew  = nan(np, length(fnames));
    grd.rh   = nan(np, length(fnames));
    grd.qv   = nan(np, length(fnames));
    grd.u    = nan(np, length(fnames));
    grd.v    = nan(np, length(fnames));
    grd.lat  = nan(np, length(fnames));
    grd.lon  = nan(np, length(fnames));
    
    % load soundings
    for ii = 1:length(fnames)
        load( fullfile( fnames(ii).folder, fnames(ii).name ) ); % snd
        wyo = getfield(load( fullfile( fnames(ii).folder, fnames(ii).name ) ), 'snd');
        %% derived quantities
        qsat = 1.0e3 * qs(snd.p*100.0, snd.t);
        
        % sonde.Specific_Humidity = 0.01*sonde.Relative_Humidity.*qsat([sonde.Temperature-273.15, sonde.Pressure]);
        if ~isfield(snd, 'qv')
            snd.qv = 0.01 * snd.rh .* qsat;
        end
        
        Pi = (snd.p*1e-3).^(Rd/Cp*(1-0.28*snd.qv*1.0e-3));
        snd.ptmp = (snd.t+273.15)./Pi;
        
        Tl = Tlcl(snd.p*1.0e2, snd.t+273.15, snd.qv);
        
        snd.etmp = snd.ptmp .* exp( (3376./Tl - 2.54) .* ...
            snd.qv*1.0e-3 .* ...
            (1 + 0.81*snd.qv*1.0e-3) ); % Bolton 1980
        
        snd.u = -snd.ws.*sind(double(snd.wd));
        snd.v = -snd.ws.*cosd(double(snd.wd));
        
        snd.thw = theta_w(snd.p*100,snd.t+273.15,snd.qv./1000);
        
        grd.time(ii) = snd.time;
        disp( datestr( snd.time(end) ) );
        
        grd.source(ii) = {'UWyo'};
        
        % average soundings to dp intervals where possible
%         pbin = [10-dp/2; dp/2+flipud(grd.p)];
        grd.temp(:,ii) = trunc(binavg( pbin, snd.p, snd.t)    );
        grd.ptmp(:,ii) = trunc(binavg( pbin, snd.p, snd.ptmp) );
        grd.etmp(:,ii) = trunc(binavg( pbin, snd.p, snd.etmp) );
        grd.thw( :,ii) = trunc(binavg( pbin, snd.p, snd.thw)  );
        grd.dew( :,ii) = trunc(binavg( pbin, snd.p, snd.td)   );
        grd.qv(  :,ii) = trunc(binavg( pbin, snd.p, snd.qv)   );
        grd.u(   :,ii) = trunc(binavg( pbin, snd.p, snd.u)    );
        grd.v(   :,ii) = trunc(binavg( pbin, snd.p, snd.v)    );
        grd.rh(  :,ii) = trunc(binavg( pbin, snd.p, snd.rh)   );
        if isfield(snd,'h')
            grd.h(  :,ii) = trunc(binavg( pbin, snd.p, snd.h)   );
        else
            grd.h( :,ii) = NaN;
        end
        if isfield(snd,'gp')
            grd.gp( :,ii) = trunc(binavg( pbin, snd.p, snd.gp)  );
        else
            grd.gp( :,ii) = NaN;
        end
        if isfield(snd,'lat')
            grd.lat(:,ii) = trunc(binavg(pbin, snd.p, snd.lat)  );
            grd.lon(:,ii) = trunc(binavg(pbin, snd.p, snd.lon)  );
        else
            grd.lat(:,ii) = NaN;
            grd.lon(:,ii) = NaN;
        end
        
        % otherwise interpolate to p grid
        [pun, un] = unique(snd.p);
        if abs(median(diff(pun))) > 2.0         
            isn = isnan(grd.temp(:,ii));
            grd.temp(isn,ii) = interp1(pun, snd.t(un)   , grd.p(isn), 'linear');
            grd.ptmp(isn,ii) = interp1(pun, snd.ptmp(un), grd.p(isn), 'linear');
            isn = isnan(grd.qv(:,ii));
            grd.etmp(isn,ii) = interp1(pun, snd.etmp(un), grd.p(isn), 'linear');
            grd.thw( isn,ii) = interp1(pun, snd.thw(un) , grd.p(isn), 'linear');
            grd.dew( isn,ii) = interp1(pun, snd.td(un)  , grd.p(isn), 'linear');
            grd.qv(  isn,ii) = interp1(pun, snd.qv(un)  , grd.p(isn), 'linear');
            isn = isnan(grd.u(:,ii));
            grd.u(   isn,ii) = interp1(pun, snd.u(un)   , grd.p(isn), 'linear');
            grd.v(   isn,ii) = interp1(pun, snd.v(un)   , grd.p(isn), 'linear');
            grd.rh(  isn,ii) = interp1(pun, snd.rh(un)  , grd.p(isn), 'linear');
            if isfield(snd,'h')
                grd.h( isn,ii) = interp1( pun, snd.h(un), grd.p(isn), 'linear');
            else
                grd.h( :,ii) = NaN;
            end
            if isfield(snd,'gp')
                grd.gp(isn,ii) = interp1( pun, snd.gp(un), grd.p(isn), 'linear');
            else
                grd.gp( :,ii) = NaN;
            end
            if isfield(snd,'lat')
                isn = isnan(grd.lat(:,ii));
                grd.lat(isn,ii) = interp1(pun, snd.lat(un), grd.p(isn), 'linear');
                grd.lon(isn,ii) = interp1(pun, snd.lon(un), grd.p(isn), 'linear');
            end
        end
        
    end
    
    save(fullfile(datadir, 'gridded', [stns_xor_wyo{istn} '_grd']), 'grd');
    
end % Wyo only stations

%% grid stations in IMD only
for istn = 1:length(stns_xor_imd)
    
    fnames = dir([matdir stns_xor_imd{istn} '_*.mat']);
    
    tmp = cell2mat({fnames.name}');
    time_string = tmp(:,end-15:end-4);
    time = datenum(time_string, 'yyyymmddHHMM');
    
    % initialize grd structure    % spd moved out of load files loop
    grd   = [];
    grd.p = (1015:-dp:10)';
    np = length(grd.p);
    grd.time = nan(1,  length(fnames));
    grd.timestr = [];
    grd.temp = nan(np, length(fnames));
    grd.ptmp = nan(np, length(fnames));
    grd.etmp = nan(np, length(fnames));
    grd.thw  = nan(np, length(fnames));
    grd.dew  = nan(np, length(fnames));
    grd.rh   = nan(np, length(fnames));
    grd.qv   = nan(np, length(fnames));
    grd.h    = nan(np, length(fnames));
    grd.gp   = nan(np, length(fnames));
    grd.u    = nan(np, length(fnames));
    grd.v    = nan(np, length(fnames));
    grd.lat  = nan(np, length(fnames));
    grd.lon  = nan(np, length(fnames));
    
    % load soundings
    for ii = 1:length(fnames)
        load( fullfile( fnames(ii).folder, fnames(ii).name ) ); % snd
        %wyo = getfield(load( fullfile( fnames(ii).folder, fnames(ii).name ) ), 'snd');
        %% derived quantities
        qsat = 1.0e3 * qs(snd.p*100.0, snd.t);
        
        % sonde.Specific_Humidity = 0.01*sonde.Relative_Humidity.*qsat([sonde.Temperature-273.15, sonde.Pressure]);
        if ~isfield(snd, 'qv')
            snd.qv = 0.01 * snd.rh .* qsat;
        end
        
        Pi = (snd.p*1e-3).^(Rd/Cp*(1-0.28*snd.qv*1.0e-3));
        snd.ptmp = (snd.t+273.15)./Pi;
        
        Tl = Tlcl(snd.p*1.0e2, snd.t+273.15, snd.qv);
        
        snd.etmp = snd.ptmp .* exp( (3376./Tl - 2.54) .* ...
            snd.qv*1.0e-3 .* ...
            (1 + 0.81*snd.qv*1.0e-3) ); % Bolton 1980
        
        snd.u = -snd.ws.*sind(double(snd.wd));
        snd.v = -snd.ws.*cosd(double(snd.wd));
        
        snd.thw = theta_w(snd.p*100,snd.t+273.15,snd.qv./1000);
%         snd.rh = 1e2*grd.qv( ip,ii) ./ nanmean(qsat);
        
        grd.time(ii) = snd.time;
        disp( datestr( snd.time(end) ) );
        
        grd.source(ii) = {snd.source};
        
        % average soundings to dp intervals where possible
%         pbin = [10-dp/2; dp/2+flipud(grd.p)];
        pbin = [10-dp/2; dp/2+flipud(grd.p)];
        grd.temp(:,ii) = trunc(binavg( pbin, snd.p, snd.t)    );
        grd.ptmp(:,ii) = trunc(binavg( pbin, snd.p, snd.ptmp) );
        grd.etmp(:,ii) = trunc(binavg( pbin, snd.p, snd.etmp) );
        grd.thw( :,ii) = trunc(binavg( pbin, snd.p, snd.thw)  );
        grd.dew( :,ii) = trunc(binavg( pbin, snd.p, snd.td)   );
        grd.qv(  :,ii) = trunc(binavg( pbin, snd.p, snd.qv)   );
        grd.u(   :,ii) = trunc(binavg( pbin, snd.p, snd.u)    );
        grd.v(   :,ii) = trunc(binavg( pbin, snd.p, snd.v)    );
        grd.rh(  :,ii) = trunc(binavg( pbin, snd.p, snd.rh)   );
        if isfield(snd, 'h')
            grd.h(  :,ii) = trunc(binavg(pbin, snd.p, snd.h)  );
        else
            grd.h(  :,ii) = NaN;
        end        
        if isfield(snd, 'gp')
            grd.gp( :,ii) = trunc(binavg(pbin, snd.p, snd.gp) );
        else
            grd.gp( :,ii) = NaN;
        end
        if isfield(snd,'lat')
            grd.lat(:,ii) = trunc(binavg(pbin, snd.p, snd.lat)  );
            grd.lon(:,ii) = trunc(binavg(pbin, snd.p, snd.lon)  );
        else
            grd.lat(:,ii) = NaN;
            grd.lon(:,ii) = NaN;
        end
        
        % otherwise interpolate to p grid
        [pun, un] = unique(snd.p);
        if abs(median(diff(pun))) > 2.0         
            isn = isnan(grd.temp(:,ii));
            grd.temp(isn,ii) = interp1(pun, snd.t(un)   , grd.p(isn), 'linear');
            grd.ptmp(isn,ii) = interp1(pun, snd.ptmp(un), grd.p(isn), 'linear');
            isn = isnan(grd.qv(:,ii));
            grd.etmp(isn,ii) = interp1(pun, snd.etmp(un), grd.p(isn), 'linear');
            grd.thw( isn,ii) = interp1(pun, snd.thw(un) , grd.p(isn), 'linear');
            grd.dew( isn,ii) = interp1(pun, snd.td(un)  , grd.p(isn), 'linear');
            grd.qv(  isn,ii) = interp1(pun, snd.qv(un)  , grd.p(isn), 'linear');
            isn = isnan(grd.u(:,ii));
            grd.u(   isn,ii) = interp1(pun, snd.u(un)   , grd.p(isn), 'linear');
            grd.v(   isn,ii) = interp1(pun, snd.v(un)   , grd.p(isn), 'linear');
            grd.rh(  isn,ii) = interp1(pun, snd.rh(un)  , grd.p(isn), 'linear');
            if isfield(snd,'h')
                grd.h( isn,ii) = interp1(binavg( pun, snd.h(un), grd.p(isn), 'linear')   );
            else
                grd.h( :,ii) = NaN;
            end
            if isfield(snd,'gp')
                grd.gp(isn,ii) = interp1(binavg( pun, snd.gp(un), grd.p(isn), 'linear')  );
            else
                grd.gp( :,ii) = NaN;
            end
            if isfield(snd,'lat')
                isn = isnan(grd.lat(:,ii));
                grd.lat(isn,ii) = interp1(pun, snd.lat(un), grd.p(isn), 'linear');
                grd.lon(isn,ii) = interp1(pun, snd.lon(un), grd.p(isn), 'linear');
            end
        end
    end
    
    save(fullfile(datadir, 'gridded', [stns_xor_imd{istn} '_grd']), 'grd');
    
end % IMD only stations