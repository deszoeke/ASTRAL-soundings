FUNCTION THETA_E, T, p, q
;
; T, Td in centigrade; p in mb
;
        r = 1000.*q
	e = q*p / (0.622 +q)
	temp = t + 273.16 ; Kelvin

;/* compute lifting condensation level temperature */

        tlcl = 2840. / (3.5 * alog( temp) - alog(e) - 4.805) + 55.

thetae = temp * (1000. / p) ^ (0.2854 *(1. - .28*.001 * r)) * $
           exp((3.376 / tlcl - 0.00254) * r * (1. + 0.81 * 0.001 * r))

return, thetae
end




