function [meanbin, stdbin, skewbin, nbin]=binavg(bin,x,X,tol)
% [meanbin, stdbin, skewbin, nbin]=binavg(bin,x,X,[tol])
% bin average 1D data on semiopen intervals x=[bin(i) bin(i+1)).
% Compute mean, standard deviation and skewness in bins.
% Returns NaNs if fewer than tol data are present. Default tol is 1. If tol<1, tol is
% treated as the fraction of the data available in the bin.
%
% Simon de Szoeke, after Jody Klymak and Sasha Perlin's bindata
% 2011-10-29

% extension possible to weighted averaging, rather than bins:
% loop over number of samples wanted in a weighted window
% compute weights based on a specified function w(x), parallel to indx definition
% multiply by weights in accumulator loop
% divide by weights in averager.

% bin in y dimension for 2d bins - make coordinates x,y a plaid grid & vectorize

% start bin averaging

if nargin<4
    tolx=1; % default tolerance
elseif tol<1 % scale fractional tolerance to number of counts
    % count how many data could be in the bin
    indxall=floor(interp1(bin,1:length(bin),x,'linear'));
    mult=histc(indxall,1:length(bin));
    tolx=tol*mult;
else % tolx specified in number of counts of data
    tolx=tol;
end

% use only finite x and X
isf=isfinite(x) & isfinite(X);
xf=x(isf);
Xf=X(isf);

indx=floor(interp1(bin,1:length(bin),xf,'linear'));
%make a special bin at the end to catch data not in any bin
indx(indx<1 | isnan(indx))=length(bin)+1;

% loop and average in bins fast
[nbin, res1,res2,res3]=deal(zeros(length(bin)+1,1));
for i=1:length(Xf) % vectorizes fast
    res1(indx(i))=res1(indx(i))+Xf(i);
    res2(indx(i))=res2(indx(i))+Xf(i)*Xf(i);
    res3(indx(i))=res3(indx(i))+Xf(i)*Xf(i)*Xf(i);

    nbin(indx(i))=nbin(indx(i))+1;
end
% averages of data not in a bin
% outbin1=res1(end);
% outbin2=res2(end);
% outbin3=res3(end);
% nout=nbin(end);
% truncate data to bins
res1=res1(1:length(bin));
res2=res2(1:length(bin));
res3=res3(1:length(bin));
nbin=nbin(1:length(bin));

% compute moments in each bin (biased estimators, always normalized by nbin)
res1(nbin<tolx)=NaN;
meanbin=res1./nbin;
mb2=meanbin.*meanbin;
varbin=res2./nbin-mb2; % 2d moment
stdbin=sqrt(varbin);
% mom3bin=(res3-3*(meanbin.*res2+mb2.*res1))./nbin - mb2.*meanbin;
mom3bin=(res3-3*(meanbin.*res2))./nbin + 2*mb2.*meanbin;
skewbin=mom3bin./(varbin.*stdbin);
