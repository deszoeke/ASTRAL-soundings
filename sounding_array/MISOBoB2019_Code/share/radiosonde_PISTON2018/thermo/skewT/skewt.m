function handles=skewt(pz,tz,rhz)
% handles=skewt(pz,tz,rhz)
% pz = pressure (hPa), tz = temperature (C), rhz = relative humidity (0-1)
%
ez=6.112.*exp(17.67.*tz./(243.5+tz));
qz=rhz.*0.622.*ez./(pz-ez);
chi=log(pz.*qz./(6.112.*(0.622+qz)));
tdz=243.5.*chi./(17.67-chi);
%
p=1050:-5:100;
pplot=transpose(p);
ps=length(p);
% t0=[-48:2:50];
t0=-23:.5:35;
ts=length(t0);
% vectorize (Simon de Szoeke)
tem=repmat(t0',[1 ps])+30*log(0.001*repmat(p,[ts 1]));
thet=(273.15+tem).*(1000./repmat(p,[ts 1])).^.287;
es=6.112.*exp(17.67.*tem./(243.5+tem));
q=622.*es./(repmat(p,[ts 1])-es);
thetaea=thet.*exp(2.5.*q./(tem+273.15));
% for i=1:ts,
%    for j=1:ps,
%       tem(i,j)=t0(i)+30.*log(0.001.*p(j));
%       thet(i,j)=(273.15+tem(i,j)).*(1000./p(j)).^.287;
%       es=6.112.*exp(17.67.*tem(i,j)./(243.5+tem(i,j)));
%       q(i,j)=622.*es./(p(j)-es);
%       thetaea(i,j)=thet(i,j).*exp(2.5.*q(i,j)./(tem(i,j)+273.15));
%    end
% end

%p=transpose(p);
t0=transpose(t0);
trang=-55:5:35;
temp=transpose(tem);
theta=transpose(thet);
thetae=transpose(thetaea);
qs=transpose(sqrt(q));
[c,h0]=contour(t0,pplot,temp,trang,'k');
hold on
set(gca,'ytick',[1000:100:100])
set(gca,'yscale','log','ydir','reverse')
set(gca,'fontweight','bold')
set(gca,'ytick',[100:100:1000])
set(gca,'ygrid','on')
[c,hth]=contour(t0,pplot,theta,24,'b');
[c,hqs]=contour(t0,pplot,qs,24,'g');
[c,hte]=contour(t0,pplot,thetae,24,'r');
set([h0(:); hth(:); hqs(:); hte(:)],'linewidth',0.25)
%tsound=30.+43.5.*log(0.001.*p);
%tsoundm=tsound-30.*log(0.001.*p);
tzm=tz-30.*log(0.001.*pz);
tdzm=tdz-30.*log(0.001.*pz);
ht=plot(tzm,pz,'k',tdzm,pz,'k--');
set(ht,'linewidth',1)
hold off
xlabel('Temperature (C)','fontweight','bold')
ylabel('Pressure (mb)','fontweight','bold')
handles=[ht; h0; hth; hqs; hte];
