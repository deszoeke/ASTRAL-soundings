function cmap = b2rcolormap(clevs)
% Sets the colormap to a high-contrast nonclashing
% blue-red colormap. clevs is the array of contour levels, or the
% scalar number of levels. (Usually, enter the length
% of the colormap + 1.)
%
% Simon de Szoeke
% 1 Dec 2005
% Thanks to Todd Mitchell at JISAO for the blue colormap.

%clevs = 10:10:100; %set your own

if length(clevs)==1
  ncols = round(clevs)-1;
else
  ncols = length(clevs)-1;
end
b = [230 255 255; 160 240 255; 80 180 255; 30 110 250; 10 50 200; 10 50 120 ]/255; 
r = fliplr(b);
n = length(b)+length(r);
cmap = interp2((1:3)',(0:n-1)/(n-1),[flipud(b); r],(1:3)',(0:ncols-1)/(ncols-1),'linear');
colormap(cmap);
