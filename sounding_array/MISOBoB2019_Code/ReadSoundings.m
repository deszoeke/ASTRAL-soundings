clear all; close all;

%% Read IMD Soundings%%
kts2ms = 0.51444444;

%% set paths
sysType = computer;
username = char(java.lang.System.getProperty('user.name'));
if strncmp(sysType,'MACI64',7)     % set Mac paths
    switch username
        case 'sdeszoek'
            datadir = '/Users/sdeszoek/Data/cruises/MISOBOB_2019/SR1911/radiosonde/';
            addpath('~/Data/cruises/MISOBOB_2019/SR1911/scientific_analysis/programs/radiosonde/');
            addpath('~/Data/cruises/MISOBOB_2019/SR1911/scientific_analysis/programs/radiosonde/radiosonde_PISTON2018/thermo/')
        otherwise
            datadir = '/Volumes/cruise/SR1911/share/data/radiosonde/';
            addpath('/Volumes/cruise/SR1911/share/scripts/radiosonde/subs/');
            addpath('/Volumes/cruise/SR1911/share/scripts/radiosonde/radiosonde_PISTON2018/thermo/')
    end
end
savedir = fullfile(datadir, 'mat/');
imd_stations={'Chennai','Kolkata','PortBlair','Thiruvananthapuram','Vishakhapatnam'};

for istn = 1:length(imd_stations)
    station=char(imd_stations(istn));
    switch station
        
        case 'Chennai'
            pth=fullfile(datadir, '/Soundings from IMD/CHENNAI/');
            fldrs=dir([pth '2019*']);
            
            for mm=1:length(fldrs)
                fname=dir(fullfile(pth, fldrs(mm).name, '/SOUNDING DATA/*.txt'));
                fprintf(1,'%s\n',fname.name)
                tmp=importdata( fullfile(fname.folder, fname.name), '\t', 3);
                
                snd=[];
                snd.time = datenum(fldrs(mm).name, 'yyyymmdd');
                snd.timestr = datestr(snd.time, 'yyyy-mm-dd HH:MM');
                snd.source = 'IMD';
                snd.p=tmp.data(:,2);
                snd.t=tmp.data(:,3);
                snd.rh=tmp.data(:,4);
                snd.td=tmp.data(:,8);
                snd.ws=tmp.data(:,5) * kts2ms; % now m/s
                snd.wd=tmp.data(:,6);
                snd.gp=tmp.data(:,7);
                
                save([savedir station '_' fname.name(1:12)], 'snd')
            end
            
        case 'Kolkata'
            pth='../../data/radiosonde/Soundings from IMD/KOLKATA/';
            fldrs=dir( fullfile(pth, '2019*') );
            
            for mm=1:length(fldrs)
                fname=dir([pth fldrs(mm).name '/SOUNDING DATA/*.txt']);
                fprintf(1,'%s\n',fname.name)
                tmp=importdata( fullfile(fname.folder, fname.name), '\t', 3);
                
                snd=[];
                snd.time = datenum(fldrs(mm).name,'yyyymmdd');
                snd.timestr = datestr(snd.time, 'yyyy-mm-dd HH:MM');
                snd.source = 'IMD';
                snd.p=tmp.data(:,1);
                snd.z=tmp.data(:,6);
                snd.t=tmp.data(:,2);
                snd.rh=tmp.data(:,3);
                snd.td=tmp.data(:,8);
                snd.ws=tmp.data(:,4);
                snd.wd=tmp.data(:,5);
                snd.gp=tmp.data(:,7);
                
                save([savedir station '_' fname.name(1:12)],'snd')
            end
            
        case {'PortBlair', 'Thiruvananthapuram', 'Vishakhapatnam'}
            switch username 
                case 'sdeszoek'
                    pth='/Users/sdeszoek/Data/cruises/MISOBOB_2019/SR1911/radiosonde/';
                otherwise
                    pth='/Volumes/cruise/SR1911/share/data/radiosonde/Soundings from IMD/';
            end
            fnames=dir([pth upper(station) '/*Summary.txt']);
            
            for mm=1:length(fnames)
                fprintf(1,'%s\n',fnames(mm).name)
                tmp=importdata([fnames(mm).folder '/' fnames(mm).name],',',4);
                
                snd=[];
                try
                    snd.time = datenum(fnames(mm).name(1:11),'yyyymmdd-HH');
                catch
                    snd.time = datenum(fnames(mm).name(15:25),'yyyymmdd-HH');
                end
                snd.timestr = datestr(snd.time, 'yyyy-mm-dd HH:MM');
                snd.source = 'IMD';
                snd.p=tmp.data(:,1);
                snd.t=tmp.data(:,2);
                snd.rh=tmp.data(:,3);
                snd.td=tmp.data(:,4);
                snd.ws=tmp.data(:,8);
                snd.wd=tmp.data(:,7);
                snd.gp=tmp.data(:,5);
                save([savedir station '_' fnames(mm).name(1:8) fnames(mm).name(10:11) '00'],'snd')
            end
            
    end % station cases
end % stations

%% Sally Ride Profiles!

pth='/Users/sdeszoek/Data/cruises/MISOBOB_2019/SR1911/radiosonde/Ride/';
flist=dir( fullfile(pth, 'raw*') );
dn = NaN(size(flist));
for mm = 1:length(flist)
    dn(mm) = datenum(flist(mm).name(17:30),'yyyymmdd_HHMM');
end
ii = dn >= datenum(2019,7,1,0,0,0);
fnames = flist(ii);

for mm=1:length(fnames)
    fprintf(1,'%s\n',fnames(mm).name)
    tmp=readtable( fullfile( fnames(mm).folder, fnames(mm).name ), 'delimiter',',:', 'headerlines',9 );
    
    snd=[];
    snd.time = datenum(fnames(mm).name(17:30),'yyyymmdd_HHMM');
    snd.timestr = datestr(snd.time, 'yyyy-mm-dd HH:MM');
    snd.source = 'SR';
    snd.h = tmp{:,2}; % height (m)
    snd.p = tmp{:,6};
    snd.t = tmp{:,7}-273.15;
    snd.ptmp = tmp{:,8};
    snd.rh=tmp{:,9};
    snd.td=tmp{:,10};
    snd.qv=tmp{:,11};
    snd.wd=tmp{:,13};
    snd.ws=tmp{:,14};
    snd.lat=tmp{:,15};
    snd.lon=tmp{:,16};
    snd.u  =tmp{:,17};
    snd.v  =tmp{:,18};

    save([savedir 'SallyRide_' fnames(mm).name(17:24) fnames(mm).name(26:27) '00'],'snd')
end

