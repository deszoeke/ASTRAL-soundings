function sonde = read_sonde(ncf)
% sonde = read_TGT_sonde(filename)
% Reads variables from one netcdf file into structure sonde.


sonde.Time                   = ncread(ncf, 'Time') ;
sonde.Height                 = ncread(ncf, 'Height') ;
sonde.Pressure               = ncread(ncf, 'Pressure') ;
sonde.Temperature            = ncread(ncf, 'Temperature') ;
sonde.Dew_Point              = ncread(ncf, 'Dew_Point') ;
sonde.Relative_Humidity      = ncread(ncf, 'Relative_Humidity') ;
sonde.Wind_Speed             = ncread(ncf, 'Wind_Speed') ;
sonde.Wind_Direction         = ncread(ncf, 'Wind_Direction') ;
sonde.Ascent_rate_of_balloon = ncread(ncf, 'Ascent_rate_of_balloon') ; 
sonde.Latitude_of_balloon    = ncread(ncf, 'Latitude_of_balloon') ;
sonde.Longitude_of_balloon   = ncread(ncf, 'Longitude_of_balloon') ;
sonde.GPS_Height             = ncread(ncf, 'GPS_Height') ;
sonde.Release_Time           = ncread(ncf, 'Release_Time') ;
sonde.Release_Latitude       = ncread(ncf, 'Release_Latitude') ;
sonde.Release_Longitude      = ncread(ncf, 'Release_Longitude') ;
sonde.Release_Altitude       = ncread(ncf, 'Release_Altitude') ;
sonde.Release_Pressure       = ncread(ncf, 'Release_Pressure') ;
sonde.Release_Temperature    = ncread(ncf, 'Release_Temperature') ;
sonde.Release_RH             = ncread(ncf, 'Release_RH') ;
sonde.Release_wind           = ncread(ncf, 'Release_wind') ;
sonde.Release_wind_direction = ncread(ncf, 'Release_wind_direction') ;

% derived quantities
Rd = 287.04;
Cp = 1005.7;
% sonde.Specific_Humidity = 0.01*sonde.Relative_Humidity.*qsat([sonde.Temperature-273.15, sonde.Pressure]);
sonde.Specific_Humidity = 0.01*sonde.Relative_Humidity .* ...
    (1.0e3 * qs(sonde.Pressure*100.0, sonde.Temperature-273.15));
sonde.Pi = (sonde.Pressure*1e-3).^(Rd/Cp*(1-0.28*sonde.Specific_Humidity*1.0e-3));
sonde.Potential_Temp = sonde.Temperature./sonde.Pi;
Tl = Tlcl(sonde.Pressure*1.0e2, sonde.Temperature, sonde.Specific_Humidity*1.0e-3);
sonde.Equiv_Pot_Temp = sonde.Potential_Temp .* exp( (3376./Tl - 2.54) .* ...
                                sonde.Specific_Humidity*1.0e-3 .* ...
                                (1 + 0.81*sonde.Specific_Humidity*1.0e-3) ); % Bolton 1980

sonde.u = -sonde.Wind_Speed.*sind(double(sonde.Wind_Direction));
sonde.v = -sonde.Wind_Speed.*cosd(double(sonde.Wind_Direction));

% eval turbulence structure on finest native vertical levels
[sonde.thsort , sonde.zsth , sonde.dth , sonde.ksth ] = Thorpe_sort(sonde.Potential_Temp, sonde.Height);
[sonde.thesort, sonde.zsthe, sonde.dthe, sonde.ksthe] = Thorpe_sort(sonde.Equiv_Pot_Temp, sonde.Height);

end
