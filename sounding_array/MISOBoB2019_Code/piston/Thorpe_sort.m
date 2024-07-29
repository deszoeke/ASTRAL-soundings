function [ts, zs, d, ks] = Thorpe_sort(t, z)
% [ts, zs, d, ks] = Thorpe_sort(t, z)
% sorts potential temperature t, calculating density sorted heights
% zs = z(ks), ts = t(ks), and overturn displacements d = z - zs

[ts, ks] = sort(t);
zs = z(ks);
d = z - zs;