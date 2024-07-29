function dudx = partiald(x,u,stndim)
% dudx = partiald(x,u,stndim) is the partial derivative
% of u(x,y) defined at locations (x,y). Calculates the dudx
% in a least squares sense. The n>=3 stations are arranged along the
% dimension stndim of x and u.
xp = x - nanmean(x, stndim);
up = u - nanmean(u, stndim);
dudx = nanmean(xp.*up, stndim)./nanvar(xp, 1, stndim);
end