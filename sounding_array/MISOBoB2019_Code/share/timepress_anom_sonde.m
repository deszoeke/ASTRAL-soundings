function [ax, cb] = timepress_anom_sonde(arr, var, titlestr, levs)
% plot timeheights and anomalies
%
% Simon de Szoeke :: MISOBOB 2019
anom = @(x) x-nanmean(x,2);

%b2rcolormap(13);
panel = 1:4;
for ist = 1:4
    arr(ist).grd.timestr = datestr(arr(ist).grd.time, 'yyyy-mm-dd HH:MM'); % fill in if it's missing

    ax(ist)=subplot(4,1,panel(ist));
    
    iimd = cellfun( @(X) strcmp(X,'IMD'), arr(ist).grd.source(:));
    pcolor( arr(ist).grd.time-datenum(2019,0,0), arr(ist).grd.p, anom(arr(ist).grd.(var)) ); shading flat;
    caxis([-3 3]); axis ij
    hold on
    contour( arr(ist).grd.time-datenum(2019,0,0), arr(ist).grd.p, arr(ist).grd.(var), levs, 'k' );
    xlim([189 floor(now)-datenum(2019,0,0)]);
    ylim([100 1000]);
    set(gca, 'ytick', 100:200:1000)
    set(gca, 'xtick', (189:(floor(now)-datenum(2019,0,0))),'tickdir','out')
    cb(ist) = colorbar();
    plot(arr(ist).grd.time(iimd)-datenum(2019,0,0), 200+zeros(sum(iimd)), 'g+')
    if ist==1
        title({titlestr, arr(ist).grd.station})
    else
        title(arr(ist).grd.station)
        set(cb(ist),'visible','off')
    end
end
end