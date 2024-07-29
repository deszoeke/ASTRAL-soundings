load p.out
load porig.out
load tdifrev.out
load tdifpseudo.out
[name,time,month,date,year]=textread('header.txt','%s %f %f %f %f');
times=cellstr(num2str(time));
months=cellstr(num2str(month));
dates=cellstr(num2str(date));
years=cellstr(num2str(year));
tit=[char(name) '   ' char(times) 'Z ' char(dates) ' ' char(months) ' ' char(years)];
clf
load modsound.txt
pr=modsound(:,1);
tc=modsound(:,2);
rh=modsound(:,3);
rh=min(rh,1);
rh=max(rh,0);
success=skewt(pr,tc,rh);
title(tit)
pause
pcolor(porig,p,tdifrev);
shading interp
set(gca,'ydir','reverse')
set(gca,'xdir','reverse')
set(gca,'fontweight','bold')
xlabel('Parcel Origin Pressure (mb)','fontweight','bold')         
ylabel('Pressure to which Parcel is Lifted (mb)','fontweight','bold')
tit1=[cellstr(char('Mean Reversible Density Temperature Difference (K)'));cellstr(tit)];
title(tit1)
hold on
[c,h]=contour(porig,p,tdifrev,'k');
clabel(c,h)
hold off
pause
clf
pcolor(porig,p,tdifpseudo);
shading interp
set(gca,'ydir','reverse')
set(gca,'xdir','reverse')
set(gca,'fontweight','bold')
xlabel('Parcel Origin Pressure (mb)','fontweight','bold')         
ylabel('Pressure to which Parcel is Lifted (mb)','fontweight','bold')
tit1=[cellstr(char('Mean Pseudo-adiabatic Density Temperature Difference (K)'));cellstr(tit)];
title(tit1)
hold on
[c,h]=contourf(porig,p,tdifpseudo,'k');
clabel(c,h)
hold off

