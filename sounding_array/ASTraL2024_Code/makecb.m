% Create nice-looking colorbar
% Alex Kinsella, November 2023

function cb = makecb(cminmax,cm,titlestr,sz,ax)


cb = colorbar;
caxis(cminmax)
cb.TickLabelInterpreter = 'latex';
cb.FontSize = sz;
ylabel(cb,titlestr,'Interpreter','latex','FontSize',sz);

if nargin == 4 
    colormap(cm)
elseif nargin > 4
    colormap(ax,cm)
end
