%% Read IMD Soundings%%
clear all; close all;

savedir='../../data/radiosonde/mat/';
imd_stations={'Chennai','Kolkata','PortBlair','Thiruvananthapuram','Vishakhapatnam'};

for istn=length(imd_stations)
    station=char(imd_stations(istn));
    switch station
        
        case 'Chennai'
            pth='../../data/radiosonde/Soundings from IMD/CHENNAI/';
            fldrs=dir([pth '2019*']);
            
            for mm=1:length(fldrs)
                fname=dir(fullfile(pth, fldrs(mm).name, '/SOUNDING DATA/*.txt'));
                fprintf(1,'%s\n',fname.name)
                tmp=importdata( fullfile(fname.folder, fname.name), '\t', 3);
                
                snd=[];
                snd.time = datenum(fldrs(mm).name, 'yyyymmdd');
                snd.timestr = datestr(snd.time, 'yyyy-mm-dd HH:MM');
                snd.p=tmp.data(:,2);
                snd.t=tmp.data(:,3);
                snd.rh=tmp.data(:,4);
                snd.td=tmp.data(:,8);
                snd.ws=tmp.data(:,5);
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
            pth='../../data/radiosonde/Soundings from IMD/';
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

pth='../../data/radiosonde/Ride/';
fnames=dir( fullfile(pth, 'raw*') );

for mm=1:length(fnames)
    fprintf(1,'%s\n',fnames(mm).name)
    tmp=importdata( fullfile( fnames(mm).folder, fnames(mm).name ), ',', 9 );
    
    snd=[];
    snd.time = datenum(fnames(mm).name(17:30),'yyyymmdd_HHMM');
    snd.timestr = datestr(snd.time, 'yyyy-mm-dd HH:MM');
    snd.p=tmp.data(:,1);
    snd.t=tmp.data(:,2)-273.15;
    snd.rh=tmp.data(:,4);
    snd.td=tmp.data(:,5)-273.15;
    snd.qv=tmp.data(:,6);
    snd.ws=tmp.data(:,9);
    snd.wd=tmp.data(:,8);
    snd.lat=tmp.data(:,10);
    snd.lon=tmp.data(:,11);
    save([savedir 'SallyRide_' fnames(mm).name(17:24) fnames(mm).name(26:27) '00'],'snd')
end

