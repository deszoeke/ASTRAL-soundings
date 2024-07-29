function esat, T

; Saturation vapor pressure as a function of temperature IN C

if( max(abs(T)) gt 170 ) then begin
	print, 'esat: T must be in centigrade ',max(abs(T))
        stop
	return, 0
endif

Tref = 273.15

;Cpv = 1850.		; water vapor heat capacity at const p
;Rv = 461.51		; Gas constant for vapor
;Lref = 2.501e6 		; vaporization heat at Tref
;Cl  =  4218.		; liquid water heat capacity
;eref = 6.106		; sat p at Tref (mb)
;
;L = Lref - (Cl - Cpv)*T		; T IN CENTIGRADE
;
; This one is from Raymond & Blyth
; esat = eref * (Tref/(T +Tref))^((Cl - Cpv)/Rv) * exp( Lref/Rv/Tref - L    /Rv/(T +Tref) )
; 


; This one is from Emanuel page 116 (4.4.13)
;

ESAT = exp( 53.67957 - 6743.769/(T+Tref) - 4.8451*alog(T+Tref) )

return, esat
end
