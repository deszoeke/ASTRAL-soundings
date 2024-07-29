%% Grid Radiosondes %%
clear all; close all;

matdir='../../data/radiosonde/mat/';
stations={'Chennai','Kolkata','PortBlair','Thiruvananthapuram','Vishakhapatnam','SallyRide'};

for m=1:length(stations)
    fnames=dir([matdir char(stations(m)) '*']);
    
    dp=5;
    grd=[];
    grd.p=[1015:-dp:10];
    grd.temp=nan(length(grd.p),length(fnames));
    grd.ptmp=nan(length(grd.p),length(fnames));
    grd.etmp=nan(length(grd.p),length(fnames));
    grd.thw=nan(length(grd.p),length(fnames));
    grd.dew=nan(length(grd.p),length(fnames));
    grd.rh=nan(length(grd.p),length(fnames));
    grd.qv=nan(length(grd.p),length(fnames));
    grd.u=nan(length(grd.p),length(fnames));
    grd.v=nan(length(grd.p),length(fnames));
    grd.lat=nan(length(grd.p),length(fnames));
    grd.lon=nan(length(grd.p),length(fnames));
    
    for mm=1:length(fnames);
        load([fnames(mm).folder '/' fnames(mm).name]);
        % derived quantities
        Rd = 287.04;
        Cp = 1005.7;
        
        qsat=1.0e3 * qs(snd.p*100.0, snd.t);
        
        % sonde.Specific_Humidity = 0.01*sonde.Relative_Humidity.*qsat([sonde.Temperature-273.15, sonde.Pressure]);
        if ~isfield(snd,'qv')
            snd.qv = 0.01*snd.rh .* qsat;
        else
            snd.qv = snd.qv;
        end
        
        Pi = (snd.p*1e-3).^(Rd/Cp*(1-0.28*snd.qv*1.0e-3));
        snd.ptmp = (snd.t+273.15)./Pi;
        
        Tl = Tlcl(snd.p*1.0e2, snd.t+273.15, snd.qv);
        
        snd.etmp = snd.ptmp .* exp( (3376./Tl - 2.54) .* ...
            snd.qv*1.0e-3 .* ...
            (1 + 0.81*snd.qv*1.0e-3) ); % Bolton 1980
        
        snd.u = -snd.ws.*sind(double(snd.wd));
        snd.v = -snd.ws.*cosd(double(snd.wd));
        
        snd.thw = theta_w(snd.p*100,snd.t+273.15,snd.qv./1000)
        
        for mmm=1:length(grd.p)
            kp=find(snd.p>=grd.p(mmm)-dp/2 & snd.p<=grd.p(mmm)+dp/2);
            grd.temp(mmm,mm)=nanmean(snd.t(kp));
            grd.ptmp(mmm,mm)=nanmean(snd.ptmp(kp));
            grd.etmp(mmm,mm)=nanmean(snd.etmp(kp));
            grd.thw(mmm,mm)=nanmean(snd.thw(kp));
            grd.dew(mmm,mm)=nanmean(snd.td(kp));
            grd.qv(mmm,mm)=nanmean(snd.qv(kp));
            grd.u(mmm,mm)=nanmean(snd.u(kp));
            grd.v(mmm,mm)=nanmean(snd.v(kp));
                grd.rh(mmm,mm)= grd.qv(mmm,mm)./nanmean(qsat(kp))./0.01;

            if isfield(snd,'lat')
                grd.lat(mmm,mm)=nanmean(snd.lat(kp));
                grd.lon(mmm,mm)=nanmean(snd.lon(kp));
            end
            
        end
    end
    
    save(['../../data/radiosonde/gridded/' char(stations(m)) '_grd'],'grd');
end

%% Grid Radiosondes %%
clear all; close all;

matdir='../../data/radiosonde/mat/wyoming/';
stations={'Chennai','Kolkata','PortBlair','Vishakhapatnam','Bhubaneswar'};

for m=1:length(stations)
    fnames=dir([matdir char(stations(m)) '*']);
    
    dp=5;
    grd=[];
    grd.p=[1015:-dp:10];
    grd.temp=nan(length(grd.p),length(fnames));
    grd.ptmp=nan(length(grd.p),length(fnames));
    grd.etmp=nan(length(grd.p),length(fnames));
    grd.thw=nan(length(grd.p),length(fnames));
    grd.dew=nan(length(grd.p),length(fnames));
    grd.rh=nan(length(grd.p),length(fnames));
    grd.qv=nan(length(grd.p),length(fnames));
    grd.u=nan(length(grd.p),length(fnames));
    grd.v=nan(length(grd.p),length(fnames));
    grd.lat=nan(length(grd.p),length(fnames));
    grd.lon=nan(length(grd.p),length(fnames));
    
    for mm=1:length(fnames);
        load([fnames(mm).folder '/' fnames(mm).name]);
        % derived quantities
        Rd = 287.04;
        Cp = 1005.7;
        
        qsat=1.0e3 * qs(snd.p*100.0, snd.t);
        
        % sonde.Specific_Humidity = 0.01*sonde.Relative_Humidity.*qsat([sonde.Temperature-273.15, sonde.Pressure]);
        if ~isfield(snd,'qv')
            snd.qv = 0.01*snd.rh .* qsat;
        else
            snd.qv = snd.qv;
        end
        
        Pi = (snd.p*1e-3).^(Rd/Cp*(1-0.28*snd.qv*1.0e-3));
        snd.ptmp = (snd.t+273.15)./Pi;
        
        Tl = Tlcl(snd.p*1.0e2, snd.t+273.15, snd.qv);
        
        snd.etmp = snd.ptmp .* exp( (3376./Tl - 2.54) .* ...
            snd.qv*1.0e-3 .* ...
            (1 + 0.81*snd.qv*1.0e-3) ); % Bolton 1980
        
        snd.u = -snd.ws.*sind(double(snd.wd));
        snd.v = -snd.ws.*cosd(double(snd.wd));
        
        snd.thw = theta_w(snd.p*100,snd.t+273.15,snd.qv./1000)
        
        for mmm=1:length(grd.p)
            kp=find(snd.p>=grd.p(mmm)-dp/2 & snd.p<=grd.p(mmm)+dp/2);
            grd.temp(mmm,mm)=nanmean(snd.t(kp));
            grd.ptmp(mmm,mm)=nanmean(snd.ptmp(kp));
            grd.etmp(mmm,mm)=nanmean(snd.etmp(kp));
            grd.thw(mmm,mm)=nanmean(snd.thw(kp));
            grd.dew(mmm,mm)=nanmean(snd.td(kp));
            grd.qv(mmm,mm)=nanmean(snd.qv(kp));
            grd.u(mmm,mm)=nanmean(snd.u(kp));
            grd.v(mmm,mm)=nanmean(snd.v(kp));
                grd.rh(mmm,mm)= grd.qv(mmm,mm)./nanmean(qsat(kp))./0.01;

            if isfield(snd,'lat')
                grd.lat(mmm,mm)=nanmean(snd.lat(kp));
                grd.lon(mmm,mm)=nanmean(snd.lon(kp));
            end
            
        end
    end
    
    save(['../../data/radiosonde/gridded/' char(stations(m)) '_grd_wyoming'],'grd');
end