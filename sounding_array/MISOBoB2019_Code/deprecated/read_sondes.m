datapath = '~/Data/cruises/MISOBOB_2019/SR1911/radiosonde/';
dfiles = dir( fullfile(datapath, 'raw_data_second_*.txt') );

curly = @(X, varargin) X{varargin};
% converts dates to matlab datenumber
dateconv = @(y) cellfun( @(x) datenum(['0-0-0 ' x]), y ); % type instability of cell bites performance
nzmax = 7200;

% get attributes from first file
i = 1;
filename = fullfile([datapath, dfiles(i).name]);
dat = readtable(filename); % much faster (than textscan) once matlab compiles readtable
% metadata
datestamp = datenum(dfiles(i).name(17:end-1), 'yyyymmdd');
Variables = dat(1,:); % a table named by variables and their units
VariableNames = Variables.Properties.VariableNames;
varj = [2 4:length(VariableNames)];

% initialize big data arrays in struct sonde
for jv = varj
    varname = VariableNames{jv};
    sonde.(varname) = NaN(nzmax, length(dfiles));
end

% load data from each file
for i = 1:length(dfiles)
    filename = fullfile([datapath, dfiles(i).name]);
    dat = readtable(filename); % much faster (than textscan) once matlab compiles readtable
%     for it = 2:length(dat.TimeUTC)
%         sonde.mtime(it,i) = datenum(['0-0-0 ' dat.TimeUTC{it}]);
%     end
    % just use the launch time
    sonde.mtime = datestamp + dateconv( dat.TimeUTC(2) );
    nz = length(dat.P) - 1;
    
    % read all the variables into sonde.(varname)(height,launch)
    % where the variables are taken from the Vaisala file
    for jv = varj
        varname = VariableNames{jv};        
        sonde.(varname)(1:nz,i) = single( cellfun(@str2num,dat.(varname)(2:end)) );
    end
end

% prototype read from Andrea Jenney
% sonde(i).pressure = single( cellfun(@str2num,dat.P(2:end)) );

% derived quantities
Rd = 287.04;
Cp = 1005.7;
% sonde.Specific_Humidity = 0.01*sonde.Relative_Humidity.*qsat([sonde.Temp-273.15, sonde.P]);
% sonde.Pi = (sonde.Pressure*1e-3).^(Rd/Cp*(1-0.28*sonde.SpecHum*1.0e-3));
sonde.Pi = sonde.Temp ./ sonde.PotTemp;
Tl = Tlcl(sonde.P*1.0e2, sonde.Temp, sonde.SpHum*1.0e-3);
sonde.EquivPotTemp = sonde.PotTemp .* exp( (3376./Tl - 2.54) .* sonde.SpHum*1.0e-3 .* ...
                                (1 + 0.81*sonde.SpHum*1.0e-3) ); % Bolton 1980
sonde.ThetaW = theta_w(sonde.P*1.0e2, sonde.Temp, sonde.SpHum*1.0e-3);
% varname = VariableNames{4}; % 4 pressure

clf
plot(sonde.PotTemp, sonde.P, 'r')
hold on 
plot(sonde.ThetaW, sonde.P, 'b')
plot((sonde.Dewp+273.15) ./ sonde.Pi, sonde.P, 'g')
axis ij; ylim([80 1010]); set(gca,'yscale','log')
xlim([290 370])
% Thorpe sorts potential temperature, calculating density sorted heights
% zs = z(ks), ts = t(ks), and overturn displacements d = z - zs
[ts, ks] = sort(sonde.PotTemp, 1);
zs = sonde.HeightMSL(ks);
d = sonde.HeightMSL - zs;
[tse, kse] = sort(sonde.EquivPotTemp, 1);
zse = sonde.HeightMSL(kse);
de = sonde.HeightMSL - zse;




