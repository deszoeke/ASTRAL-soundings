% Compare integrated vapor WVP from microwave radiometer with soundings.

thermo_constants; % thermo functions such as qs require global constants definitions
read_pbin_all_TGT_sondes;
dp = repmat(diff(pint),[1 size(shum,2)]);
dp(end,:) = psfc-pint(end-1);
Mvap = nansum(1.e2*dp.*shum/1.e3)/g; % kg/m^2
rhol = 1.e3; % kg m^-3
wvp = Mvap / rhol * 1.e3; % mm = kg/m^2

colmfile = fullfile('/Users/sdeszoek/Data/cruises/PISTON_MISOBOB_2018/TGT/Radiometer/operational/lv2/ColInt_Obs.txt');
fid = fopen(colmfile,'r');
frewind(fid);
tmp = textscan(fid,'%f %f/%f/%f %f:%f:%f %f  %f  %f  %f %*[^\n]',...
    'multipledelimsasone',1,'delimiter',{' ',','},'headerlines',1); % 5G file takes several minutes
fclose(fid);
c=cell2mat(tmp);% Record, Date/Time, Record_Type, Int_Vapor_cm, Int_Liquid_mm, Cloud_Base_km
%n bb mm/dd/yy hh:mm:ss
year = c(:,4) + 2000;
colm.yday = datenum([0*year c(:,[2:3 5:7])]);
colm.vapr = c(:, 9); % cm
colm.liqd = c(:,10); % mm = kg/m^2
colm.base = c(:,11);

semilogy(colm.yday,colm.vapr*10)
hold on
plot(colm.yday,colm.liqd)
ylabel('kg m^{-2}')
xlabel('yearday 2018')
legend('vapor','liquid')
title('PISTON 2018, microwave radiometer water paths (junk)','fontweight','normal')
set(gca,'fontsize',16)
saveas(gcf,'PISTON2018_colintMWR.png')
% MWR useless without reprocessing.